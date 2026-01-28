# ============================================================
# Figure 2. Adjacency matrices from estimated MMI networks for
# MACE and non-MACE groups.
#
# Panels A and B: GGM-based networks for MACE and non-MACE.
# Panels C and D: GGMridge-based networks.
# Panels E and F: EBICglasso-based networks.
# ============================================================

# ---- Libraries (keep as you used; minimal changes) ----
library(GeneNet)
library(corrplot)
library(dplyr)
library(GGMridge)
library(huge)

# ---- File paths (edit to match your repo structure) ----
data_path <- "data/NonPrevExtMaceSamples_IDedDedupMtbs.csv"
ext_mace_ids_path <- "data/ExtMacePatids.txt"
no_ext_mace_ids_path <- "data/NoExtMacePatids.txt"

# ---- Output folder (optional, but recommended) ----
out_dir <- "outputs/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read data + patient ID lists ----
data <- read.csv(data_path)
ext_mace_patids <- readLines(ext_mace_ids_path)
no_ext_mace_patids <- readLines(no_ext_mace_ids_path)

# ============================================================
# Split into groups
# ============================================================
ext_mace_subset <- data %>% filter(patid %in% ext_mace_patids)
no_ext_mace_subset <- data %>% filter(patid %in% no_ext_mace_patids)

# Use dplyr::select to avoid select() conflicts
ext_mace_metabolites_data <- ext_mace_subset %>% dplyr::select(-patid)
no_ext_mace_metabolites_data <- no_ext_mace_subset %>% dplyr::select(-patid)

p <- ncol(ext_mace_metabolites_data)  # should be 614

# ============================================================
# Panel A/B: GGM (GeneNet) adjacency via edge testing
# ============================================================

ext_mace_pcor_matrix <- ggm.estimate.pcor(ext_mace_metabolites_data)
Pcor_AllEdges_MACEgroup <- network.test.edges(ext_mace_pcor_matrix, fdr = TRUE, plot = FALSE)

no_ext_mace_pcor_matrix <- ggm.estimate.pcor(no_ext_mace_metabolites_data)
Pcor_AllEdges_noMACEgroup <- network.test.edges(no_ext_mace_pcor_matrix, fdr = TRUE, plot = FALSE)

# Keep high-confidence edges (as in your code)
Pcor_AllEdges_MACEgroup_lim <- Pcor_AllEdges_MACEgroup %>% filter(prob > 0.9)
Pcor_AllEdges_noMACEgroup_lim <- Pcor_AllEdges_noMACEgroup %>% filter(prob > 0.9)

parcormatrix_mace <- matrix(0, nrow = p, ncol = p)
if (nrow(Pcor_AllEdges_MACEgroup_lim) > 0) {
  for (i in 1:nrow(Pcor_AllEdges_MACEgroup_lim)) {
    row <- Pcor_AllEdges_MACEgroup_lim[i, 2]
    column <- Pcor_AllEdges_MACEgroup_lim[i, 3]
    parcormatrix_mace[row, column] <- 1
    parcormatrix_mace[column, row] <- 1
  }
}

parcormatrix_nomace <- matrix(0, nrow = p, ncol = p)
if (nrow(Pcor_AllEdges_noMACEgroup_lim) > 0) {
  for (i in 1:nrow(Pcor_AllEdges_noMACEgroup_lim)) {
    row <- Pcor_AllEdges_noMACEgroup_lim[i, 2]
    column <- Pcor_AllEdges_noMACEgroup_lim[i, 3]
    parcormatrix_nomace[row, column] <- 1
    parcormatrix_nomace[column, row] <- 1
  }
}

# Save Panel A (KEEP original corrplot settings)
pdf(file.path(out_dir, "Fig2A_GGM_MACE.pdf"), width = 7, height = 7)
corrplot(parcormatrix_mace, method = "color", tl.col = "white", tl.cex = 0.001)
dev.off()

# Save Panel B (KEEP original corrplot settings)
pdf(file.path(out_dir, "Fig2B_GGM_noMACE.pdf"), width = 7, height = 7)
corrplot(parcormatrix_nomace, method = "color", tl.col = "white", tl.cex = 0.001)
dev.off()

# ============================================================
# Panel C/D: GGMridge-based adjacency
# (Your original PCMRidge() retained, only wrapped for saving)
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
  }
  for (i in 1:pdat) {
    lambda_numerator <- lambda_numerator + var_estimate(i, i)
  }
  
  lambda_denominator <- 0
  for (i in 1:pdat) {
    for (j in 1:pdat) {
      if (i != j) lambda_denominator <- lambda_denominator + S[i, j]^2
    }
  }
  for (i in 1:pdat) {
    lambda_denominator <- lambda_denominator + (S[i, i] - 1)^2
  }
  
  lambda_star <- lambda_numerator / lambda_denominator
  lambda_0 <- lambda_star / (1 - lambda_star)
  
  lambda.array <- seq(from = lambda_0, to = lambda_0 + 3, length = 10)
  pcut.array <- seq(from = 0.01/(pdat^2), to = 0.05, length = 10)
  
  tpe <- lambda.pcut.cv(x = dat, lambda = lambda.array, pcut = pcut.array, fold = 10)
  w.mintpe <- which(tpe == min(tpe), arr.ind = TRUE)
  lambda <- lambda.array[w.mintpe[1L]]
  alpha <- pcut.array[w.mintpe[2L]]
  
  partial <- solve(S + lambda * diag(pdat))
  partial <- -scaledMat(x = partial)
  
  w.upper <- which(upper.tri(diag(pdat)))
  psi <- transFisher(x = partial[w.upper])
  
  efron.fit <- getEfronp(z = psi)
  psi_standardized <- (psi - efron.fit$mu0hat) / (efron.fit$sigma0hat)
  p_value <- 2 * (1 - pnorm(abs(psi_standardized)))
  
  edge3 <- ifelse(p_value < alpha, 1, 0)
  
  edge_matrix <- matrix(0, nrow = pdat, ncol = pdat)
  edge_matrix[upper.tri(edge_matrix)] <- edge3
  edge_matrix <- edge_matrix + t(edge_matrix)
  edge_matrix
}

parcormatrix_mace <- PCMRidge(ext_mace_metabolites_data)
parcormatrix_nomace <- PCMRidge(no_ext_mace_metabolites_data)

# Save Panel C (KEEP original corrplot settings)
pdf(file.path(out_dir, "Fig2C_GGMridge_MACE.pdf"), width = 7, height = 7)
corrplot(parcormatrix_mace, method = "color", tl.col = "white", tl.cex = 0.001)
dev.off()

# Save Panel D (KEEP original corrplot settings)
pdf(file.path(out_dir, "Fig2D_GGMridge_noMACE.pdf"), width = 7, height = 7)
corrplot(parcormatrix_nomace, method = "color", tl.col = "white", tl.cex = 0.001)
dev.off()

# ============================================================
# Panel E/F: EBICglasso (huge)
# ============================================================
  
# MACE
x_mace <- as.matrix(ext_mace_metabolites_data)
fit_mace <- huge(x_mace, method = "glasso", nlambda = 30)
select_ebic_mace <- huge.select(fit_mace, criterion = "ebic")
parcormatrix_mace <- select_ebic_mace$refit
if (length(dim(parcormatrix_mace)) == 3) parcormatrix_mace <- parcormatrix_mace[, , 1]

# non-MACE
x_nomace <- as.matrix(no_ext_mace_metabolites_data)
fit_nomace <- huge(x_nomace, method = "glasso", nlambda = 30)
select_ebic_nomace <- huge.select(fit_nomace, criterion = "ebic")
parcormatrix_nomace <- select_ebic_nomace$refit
if (length(dim(parcormatrix_nomace)) == 3) parcormatrix_nomace <- parcormatrix_nomace[, , 1]

# Save Panel E (KEEP original corrplot settings)
pdf(file.path(out_dir, "Fig2E_EBICglasso_MACE.pdf"), width = 7, height = 7)
corrplot(parcormatrix_mace, method = "color", tl.col = "white", tl.cex = 0.001)
dev.off()

# Save Panel F (KEEP original corrplot settings)
pdf(file.path(out_dir, "Fig2F_EBICglasso_noMACE.pdf"), width = 7, height = 7)
corrplot(parcormatrix_nomace, method = "color", tl.col = "white", tl.cex = 0.001)
dev.off()
