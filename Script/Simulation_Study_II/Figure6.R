
# 1. Set up----

#install.packages("huge")
library(huge)
#install.packages("semTools")
library(semTools)
#install.packages("corrplot")
library(corrplot)
#install.packages("dplyr")
library(dplyr)
#install.packages("readr")
library(readr)
#install.packages("GGMridge")
library(GGMridge)
#install.packages("GeneNet")
library(GeneNet)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("moments")
library(moments)

Ridgeseed <- c(21169, 24918, 30401, 31829, 44683, 46825, 59075, 83133, 93595, 145952,
               153263, 211035, 223212, 236829, 264699, 277048, 295482, 309734, 314778, 341850,
               343070, 361673, 386319, 418871, 434826, 456011, 457665, 463159, 516642, 517898,
               537280, 598693, 634619, 645741, 655370, 685033, 728341, 734175, 738175, 769225,
               773687, 790298, 802683, 870961, 875458, 890909, 891868, 893344, 895078, 935644)

seed <- Ridgeseed[1]
set.seed(seed)

n <- 1000  # Number of samples
p <- 100     # Number of variables
mu <- rep(0, p)
num_simulations <- 500
subsample_sizes <- c(200, 400, 600, 800, 1000)

metrics <- read.csv("2-2_metrics_values.csv")
metrics <- subset(metrics, abs(skewness) <= 5)

data_forsigma <- huge.generator(n, d = p, graph = "cluster", vis = FALSE)
Sigma <- as.matrix(data_forsigma$sigma)
diag(Sigma) <- 1
actual <- as.matrix(data_forsigma$theta)

metrics_sample <- metrics[sample(nrow(metrics), p),]

X <- mvrnonnorm(n, mu, Sigma, 
                skewness = metrics_sample$skewness, 
                kurtosis = metrics_sample$kurtosis, 
                empirical = FALSE)

# 2. Ridge Method Function----

# Partial Correlation Matrix Estimation using Ridge Penalty

PCMRidge <- function(dat) {
  ndat <- nrow(dat)
  pdat <- ncol(dat)
  
  colmeans <- colMeans(dat)
  S <- cov(dat)
  W <- S * (ndat - 1) / ndat
  
  Wkij <- array(0, dim = c(pdat, pdat, ndat))
  for (k in 1:ndat) {
    x <- as.numeric(dat[k, ] - colmeans)
    #Wkij[,,k] <- (dat[k,] - colmeans) %*% t(dat[k,] - colmeans)
    Wkij[,,k] <- x %*% t(x)
  }
  
  Wij <- apply(Wkij, 1:2, mean)
  
  # Function to compute Var_hat(s_ij) for given (i, j)
  var_estimate <- function(i, j) {
    sum_k <- sum((Wkij[i, j, ] - Wij[i, j])^2)
    result <- (ndat / (ndat - 1)^3) * sum_k
    return(result)
  }
  
  # Compute the numerator of λ⋆
  lambda_numerator <- 0
  for (i in 1:pdat) {
    for (j in 1:pdat) {
      if (i != j) {
        lambda_numerator <- lambda_numerator + var_estimate(i, j)
      }
    }
  }
  for (i in 1:pdat) {
    lambda_numerator <- lambda_numerator + var_estimate(i, i)
  }
  
  # Compute the denominator of λ⋆
  lambda_denominator <- 0
  for (i in 1:pdat) {
    for (j in 1:pdat) {
      if (i != j) {
        lambda_denominator <- lambda_denominator + S[i, j]^2
      }
    }
  }
  for (i in 1:pdat) {
    lambda_denominator <- lambda_denominator + (S[i, i] - 1)^2
  }
  
  # Calculate λ⋆
  lambda_star <- lambda_numerator / lambda_denominator
  lambda_star #0.04489079
  
  lambda_0 = lambda_star/(1-lambda_star)
  
  # Selection of a lambda and a p-value cutoff
  lambda.array <- seq(from = lambda_0, to = lambda_0+3, length = 10)
  pcut.array <- seq(from = 0.01/(p^2), to = 0.05, length = 10)
  tpe <- lambda.pcut.cv(x = dat,
                        lambda = lambda.array,
                        pcut = pcut.array,
                        fold = 10)
  w.mintpe <- which(tpe == min(tpe), arr.ind = TRUE)
  (lambda <- lambda.array[w.mintpe[1L]])
  (alpha <- pcut.array[w.mintpe[2L]])
  
  ## partial correlation matrix use ridge inverse
  partial <- solve(S + lambda * diag(pdat))
  partial <- -scaledMat(x = partial)
  
  # Fisher's Z transformation of upper diagonal of the partial correlation matrix
  w.upper <- which(upper.tri(diag(pdat)))
  psi <- transFisher(x = partial[w.upper])
  
  # get p-values from empirical null distribution
  efron.fit <- getEfronp(z = psi)
  psi_standardized <- (psi - efron.fit$mu0hat) / (efron.fit$sigma0hat)
  p_value <- 2 * (1 - pnorm(abs(psi_standardized)))
  #p_value <- 2 * (1 - pnorm(abs(psi), mean = efron.fit$mu0hat, sd = efron.fit$sigma0hat))
  edge3 <- ifelse(p_value < alpha, 1, 0)
  
  edge_matrix <- matrix(0, nrow = p, ncol = p)
  edge_matrix[upper.tri(edge_matrix)] <- edge3
  edge_matrix <- edge_matrix + t(edge_matrix)
  return(edge_matrix)
}

# 3. Run Simulations----

# Simulate 500 times
for (subsample in subsample_sizes) {
  FN_counts <- matrix(0, p, p) # Reset for each subsample size
  FP_counts <- matrix(0, p, p) # Reset for each subsample size
  
  sens <- rep(0, num_simulations)
  spec <- rep(0, num_simulations)
  PPV <- rep(0, num_simulations)
  NPV <- rep(0, num_simulations)
  
  for (i in 1:num_simulations) {
    set.seed(i)
    dat <- X[sample(nrow(X), subsample), ]
    parcormatrix <- PCMRidge(dat)
    predicted <- parcormatrix
    
    # Identify FN edges
    FN_mask <- (predicted == 0) & (actual == 1)
    FP_mask <- (predicted == 1) & (actual == 0)
    # Count occurrences of each FN edge
    FN_counts[FN_mask] <- FN_counts[FN_mask] + 1
    FP_counts[FP_mask] <- FP_counts[FP_mask] + 1
    
    # Identify FN edges
    FN <- sum((predicted == 0) & (actual == 1))/2
    FP <- sum((predicted == 1) & (actual == 0))/2
    TP <- sum((predicted == 1) & (actual == 1))/2
    TN <- sum((predicted == 0) & (actual == 0))/2
    
    sens[i] <- TP/(TP+FN)
    spec[i] <- TN/(TN+FP)
    PPV[i] <- TP/(TP+FP)
    NPV[i] <- TN/(TN+FN)
    
    print(paste("Subsample:", subsample, "Iteration:", i))
  }
  
  # Convert FN count matrix to a list of edges
  upper_tri_indices <- which(upper.tri(FN_counts), arr.ind = TRUE)
  count_data <- data.frame(
    row = upper_tri_indices[, 1],
    col = upper_tri_indices[, 2],
    FN_count = FN_counts[upper_tri_indices],
    FP_count = FP_counts[upper_tri_indices]
  )
  
  # Compute skewness and kurtosis for each variable (column)
  skewness_values <- apply(X, 2, skewness, na.rm = TRUE)
  kurtosis_values <- apply(X, 2, moments::kurtosis, na.rm = TRUE)
  mean_values <- apply(X, 2, mean, na.rm = TRUE)
  variance_values <- apply(X, 2, var, na.rm = TRUE)
  
  count_data <- count_data %>%
    mutate(
      skew1 = skewness_values[row],
      skew2 = skewness_values[col],
      kurt1 = kurtosis_values[row],
      kurt2 = kurtosis_values[col],
      mean1 = mean_values[row],
      mean2 = mean_values[col],
      var1 = variance_values[row],
      var2 = variance_values[col],
      skewavg = 0.5 * (skew1 + skew2),
      kurtavg = 0.5 * (kurt1 + kurt2)
    )
  
  metric_data <- data.frame(
    sens = sens,
    spec = spec,
    PPV = PPV,
    NPV = NPV
  )
  
  write.csv(count_data, paste0("Synthetic_count_data_", subsample,"_Ridge",seed,".csv"))
  write.csv(metric_data, paste0("Synthetic_metric_data_", subsample,"_Ridge_",seed,".csv"))
}
