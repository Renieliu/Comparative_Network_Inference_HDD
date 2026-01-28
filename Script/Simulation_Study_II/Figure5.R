
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


rm(list = ls())

# random number
GGMseed <- c(21169, 24918, 30401, 31829, 44683, 46825, 59075, 83133, 93595, 145952,
             153263, 211035, 223212, 236829, 264699, 277048, 295482, 309734, 314778, 341850,
             343070, 361673, 386319, 418871, 434826, 456011, 457665, 463159, 516642, 517898,
             537280, 598693, 634619, 645741, 655370, 685033, 728341, 734175, 738175, 769225,
             773687, 790298, 802683, 870961, 875458, 890909, 891868, 893344, 895078, 935644)

seed <- GGMseed[1] #repeat for all seed in GGMseed
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

# Simulate 500 times
for (subsample in subsample_sizes) {
  FN_counts <- matrix(0, p, p) # Reset for each subsample size
  FP_counts <- matrix(0, p, p) # Reset for each subsample size
  
  sens <- rep(0, num_simulations)
  spec <- rep(0, num_simulations)
  PPV <- rep(0, num_simulations)
  NPV <- rep(0, num_simulations)

  # Simulate 500 times
  for (i in 1:num_simulations) {
    set.seed(i)
    dat <- X[sample(nrow(X), subsample), ]
    
    pcor_matrix <- ggm.estimate.pcor(dat)
    Pcor_AllEdges <- network.test.edges(pcor_matrix, fdr = TRUE, plot = FALSE)
    Pcor_AllEdges_lim <- filter(Pcor_AllEdges, prob > 0.9)
    parcormatrix <- matrix(0, p, p)
    
    for (j in seq_len(nrow(Pcor_AllEdges_lim))) {
      row <- Pcor_AllEdges_lim[j, 2]
      column <- Pcor_AllEdges_lim[j, 3]
      parcormatrix[row, column] <- 1
      parcormatrix[column, row] <- 1
    }
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
  
  write.csv(count_data, paste0("Synthetic_count_data_", subsample,"_GGM_",seed,".csv"))
  write.csv(metric_data, paste0("Synthetic_metric_data_", subsample,"_GGM_",seed,".csv"))
}

