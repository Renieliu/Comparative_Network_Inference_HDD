# ============================================================
# Fig 4: Adjacency matrices from estimated zoom-in (50 × 50) MMI networks for synthetic multivariate normal data 
# Method: GGM (GeneNet partial correlation) + edge test (FDR) + prob threshold
#
# Output: heatmap-like adjacency matrix plot via corrplot
# ============================================================

library(corrplot)
library(dplyr)
library(readr)
library(GeneNet)
library(igraph)

subsamplesize <- 300 
seed <- 123

# ---- Create synthetic multivariate normal data ----
set.seed(seed)
pcor <- ggm.simulate.pcor(614, etaA=0.01)
sim.data <-  ggm.simulate.data(2060, pcor)

ext_mace_metabolites_data <- sim.data[1:263,]
no_ext_mace_metabolites_data <- sim.data[264:2060,]
no_ext_mace <- no_ext_mace_metabolites_data

# ---- Random subsample ----
selected_rows <- sample(nrow(no_ext_mace), subsamplesize)
no_ext_mace_subset <- no_ext_mace[selected_rows, ]
no_ext_mace_metabolites_data <- no_ext_mace_subset %>% dplyr::select(-patid)

p <- ncol(no_ext_mace_metabolites_data)  # should be 614

# ---- GGM partial correlation + edge testing ----
no_ext_mace_pcor_matrix <- ggm.estimate.pcor(no_ext_mace_metabolites_data)
Pcor_AllEdges_noMACEgroup <- network.test.edges(no_ext_mace_pcor_matrix, fdr = TRUE)
Pcor_AllEdges_noMACEgroup_lim <- Pcor_AllEdges_noMACEgroup %>% filter(prob > 0.9)

# ---- Build adjacency matrix (0/1) ----
parcormatrix_nomace <- matrix(rep(0, p * p), nrow = p, ncol = p)
for (i in 1:nrow(Pcor_AllEdges_noMACEgroup_lim)) {
  row <- Pcor_AllEdges_noMACEgroup_lim[i, 2]
  column <- Pcor_AllEdges_noMACEgroup_lim[i, 3]
  parcormatrix_nomace[row, column] <- 1
  parcormatrix_nomace[column, row] <- 1
}

# ---- Plot (and optionally save) ----
out_file <- file.path(out_dir, paste0("Fig4_synMVN_GGM_adjacency_n", subsamplesize, ".png"))
png(out_file, width = 2200, height = 2200, res = 300)
corrplot(parcormatrix_nomace, method = "color", tl.col = "white", tl.cex = 0.001)
dev.off()