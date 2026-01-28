# ============================================================
# Table 2. Comparison of node-level and network-level metrics
# between MACE and non-MACE groups across three estimation methods.
# Methods: (1) GGM (GeneNet), (2) GGMridge, (3) EBICglasso (huge)
# Output: outputs/tables/table2.csv
# ============================================================

library(dplyr)
library(readr)
library(igraph)
library(GeneNet)
library(GGMridge)
library(huge)

# ---- Paths (set to your repo structure) ----
data_path <- "data/NonPrevExtMaceSamples_IDedDedupMtbs.csv"
ext_mace_ids_path <- "data/ExtMacePatids.txt"
no_ext_mace_ids_path <- "data/NoExtMacePatids.txt"
out_dir <- "outputs/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read data ----
dat <- readr::read_csv(data_path, show_col_types = FALSE)
ext_mace_patids <- readLines(ext_mace_ids_path)
no_ext_mace_patids <- readLines(no_ext_mace_ids_path)

ext_dat <- dat %>% filter(patid %in% ext_mace_patids) %>% dplyr::select(-patid)
no_dat  <- dat %>% filter(patid %in% no_ext_mace_patids) %>% dplyr::select(-patid)

p <- ncol(ext_dat) # 614

# ============================================================
# Distances between two adjacency matrices
# ============================================================

edit_dist <- function(A1, A2) {
  stopifnot(all(dim(A1) == dim(A2)))
  sum(abs(A1 - A2)) / 2  # divide by 2 for undirected symmetric adjacency
}

adj_spec_dist <- function(A1, A2) {
  stopifnot(all(dim(A1) == dim(A2)))
  ev1 <- sort(eigen(A1, symmetric = TRUE, only.values = TRUE)$values)
  ev2 <- sort(eigen(A2, symmetric = TRUE, only.values = TRUE)$values)
  sqrt(sum((ev1 - ev2)^2))
}

lap_spec_dist <- function(A1, A2) {
  stopifnot(all(dim(A1) == dim(A2)))
  L1 <- diag(rowSums(A1)) - A1
  L2 <- diag(rowSums(A2)) - A2
  ev1 <- sort(eigen(L1, symmetric = TRUE, only.values = TRUE)$values)
  ev2 <- sort(eigen(L2, symmetric = TRUE, only.values = TRUE)$values)
  sqrt(sum((ev1 - ev2)^2))
}

normlap_spec_dist <- function(A1, A2) {
  stopifnot(all(dim(A1) == dim(A2)))
  D1 <- rowSums(A1); D2 <- rowSums(A2)
  
  # avoid divide-by-zero for isolated nodes
  invsqrtD1 <- ifelse(D1 > 0, 1 / sqrt(D1), 0)
  invsqrtD2 <- ifelse(D2 > 0, 1 / sqrt(D2), 0)
  
  L1 <- diag(D1) - A1
  L2 <- diag(D2) - A2
  NL1 <- diag(invsqrtD1) %*% L1 %*% diag(invsqrtD1)
  NL2 <- diag(invsqrtD2) %*% L2 %*% diag(invsqrtD2)
  
  ev1 <- sort(eigen(NL1, symmetric = TRUE, only.values = TRUE)$values)
  ev2 <- sort(eigen(NL2, symmetric = TRUE, only.values = TRUE)$values)
  sqrt(sum((ev1 - ev2)^2))
}

# ============================================================
# Graph metrics from adjacency
# ============================================================

calc_metrics <- function(A1, A2) {
  g1 <- igraph::graph_from_adjacency_matrix(A1, mode = "undirected", diag = FALSE)
  g2 <- igraph::graph_from_adjacency_matrix(A2, mode = "undirected", diag = FALSE)
  
  out <- c(
    close_mace      = mean(count_triangles(g1)),
    close_nomace    = mean(count_triangles(g2)),
    closeopen_mace  = mean(choose(degree(g1), 2)),
    closeopen_nomace= mean(choose(degree(g2), 2)),
    open_mace       = mean(choose(degree(g1), 2) - count_triangles(g1)),
    open_nomace     = mean(choose(degree(g2), 2) - count_triangles(g2)),
    degree_cen_mace = mean(degree(g1), na.rm = TRUE),
    degree_cen_nomace=mean(degree(g2), na.rm = TRUE),
    bw_cen_mace     = mean(betweenness(g1), na.rm = TRUE),
    bw_cen_nomace   = mean(betweenness(g2), na.rm = TRUE),
    close_cen_mace  = mean(closeness(g1), na.rm = TRUE),
    close_cen_nomace= mean(closeness(g2), na.rm = TRUE),
    edge_number_mace= sum(A1),
    edge_number_nomace=sum(A2),
    ED              = edit_dist(A1, A2),
    ASD             = adj_spec_dist(A1, A2),
    LSD             = lap_spec_dist(A1, A2),
    NLSD            = normlap_spec_dist(A1, A2)
  )
  out
}

# ============================================================
# 1) GGM adjacency (GeneNet)
# ============================================================

adj_ggm <- function(x, prob_thr = 0.9) {
  pcor <- GeneNet::ggm.estimate.pcor(x)
  edges <- GeneNet::network.test.edges(pcor, fdr = TRUE, plot = FALSE)
  edges <- edges %>% filter(prob > prob_thr)
  
  A <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  if (nrow(edges) > 0) {
    for (i in seq_len(nrow(edges))) {
      r <- edges[i, 2]; c <- edges[i, 3]
      A[r, c] <- 1; A[c, r] <- 1
    }
  }
  A
}

# ============================================================
# 2) GGMridge adjacency
# ============================================================

PCMRidge <- function(dat) {
  ndat <- nrow(dat)
  pdat <- ncol(dat)
  
  colmeans <- colMeans(dat)
  S <- cov(dat)
  
  Wkij <- array(0, dim = c(pdat, pdat, ndat))
  for (k in 1:ndat) {
    x <- as.numeric(dat[k, ] - colmeans)
    Wkij[, , k] <- x %*% t(x)
  }
  Wij <- apply(Wkij, 1:2, mean)
  
  var_estimate <- function(i, j) {
    sum_k <- sum((Wkij[i, j, ] - Wij[i, j])^2)
    (ndat / (ndat - 1)^3) * sum_k
  }
  
  lambda_numerator <- 0
  for (i in 1:pdat) {
    for (j in 1:pdat) {
      if (i != j) lambda_numerator <- lambda_numerator + var_estimate(i, j)
    }
    lambda_numerator <- lambda_numerator + var_estimate(i, i)
  }
  
  lambda_denominator <- 0
  for (i in 1:pdat) {
    for (j in 1:pdat) {
      if (i != j) lambda_denominator <- lambda_denominator + S[i, j]^2
    }
    lambda_denominator <- lambda_denominator + (S[i, i] - 1)^2
  }
  
  lambda_star <- lambda_numerator / lambda_denominator
  lambda_0 <- lambda_star / (1 - lambda_star)
  
  lambda.array <- seq(from = lambda_0, to = lambda_0 + 3, length = 10)
  pcut.array <- seq(from = 0.01/(pdat^2), to = 0.05, length = 10)
  
  tpe <- GGMridge::lambda.pcut.cv(x = dat, lambda = lambda.array, pcut = pcut.array, fold = 10)
  w.mintpe <- which(tpe == min(tpe), arr.ind = TRUE)
  lambda <- lambda.array[w.mintpe[1L]]
  alpha <- pcut.array[w.mintpe[2L]]
  
  partial <- solve(S + lambda * diag(pdat))
  partial <- -GGMridge::scaledMat(x = partial)
  
  w.upper <- which(upper.tri(diag(pdat)))
  psi <- GGMridge::transFisher(x = partial[w.upper])
  
  efron.fit <- GGMridge::getEfronp(z = psi)
  psi_standardized <- (psi - efron.fit$mu0hat) / (efron.fit$sigma0hat)
  p_value <- 2 * (1 - pnorm(abs(psi_standardized)))
  
  edge3 <- ifelse(p_value < alpha, 1, 0)
  A <- matrix(0, nrow = pdat, ncol = pdat)
  A[upper.tri(A)] <- edge3
  A <- A + t(A)
  A
}

# ============================================================
# 3) EBICglasso adjacency (huge)
# ============================================================

adj_ebicglasso <- function(x, nlambda = 30) {
  fit <- huge::huge(as.matrix(x), method = "glasso", nlambda = nlambda)
  sel <- huge::huge.select(fit, criterion = "ebic")
  A <- sel$refit
  if (length(dim(A)) == 3) A <- A[, , 1]
  A
}

# ============================================================
# Bootstrap wrapper for one method
# ============================================================

bootstrap_table2 <- function(method_name, make_adj, X_mace, X_nomace, B = 100, seed = 1) {
  set.seed(seed)
  
  # point estimate
  A1 <- make_adj(X_mace)
  A2 <- make_adj(X_nomace)
  point <- calc_metrics(A1, A2)
  
  # bootstrap
  n1 <- nrow(X_mace); n2 <- nrow(X_nomace)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(point))
  colnames(boot_mat) <- names(point)
  
  for (b in seq_len(B)) {
    samp1 <- X_mace[sample.int(n1, n1, replace = TRUE), , drop = FALSE]
    samp2 <- X_nomace[sample.int(n2, n2, replace = TRUE), , drop = FALSE]
    A1b <- make_adj(samp1)
    A2b <- make_adj(samp2)
    boot_mat[b, ] <- calc_metrics(A1b, A2b)
    message(method_name, " bootstrap: ", b, "/", B)
  }
  
  ci_lo <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_hi <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  out <- data.frame(
    Method = method_name,
    Metric = names(point),
    Estimate = as.numeric(point),
    CI_low = as.numeric(ci_lo),
    CI_high = as.numeric(ci_hi),
    row.names = NULL
  )
  out
}

# ============================================================
# Run all three methods
# ============================================================

res_ggm <- bootstrap_table2(
  method_name = "GGM",
  make_adj = function(x) adj_ggm(x, prob_thr = 0.9),
  X_mace = ext_dat, X_nomace = no_dat,
  B = 100, seed = 1
)

res_ridge <- bootstrap_table2(
  method_name = "GGMridge",
  make_adj = PCMRidge,
  X_mace = ext_dat, X_nomace = no_dat,
  B = 100, seed = 1
)

res_ebic <- bootstrap_table2(
  method_name = "EBICglasso",
  make_adj = function(x) adj_ebicglasso(x, nlambda = 30),
  X_mace = ext_dat, X_nomace = no_dat,
  B = 100, seed = 1
)

table2_long <- bind_rows(res_ggm, res_ridge, res_ebic)

# Save long-format (clean + easy to reshape in LaTeX)
readr::write_csv(table2_long, file.path(out_dir, "table2_long.csv"))
