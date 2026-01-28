# ============================================================
# Figure 1: Skewness vs. Kurtosis of metabolite features
# Data: Diabetes metabolomics dataset (614 metabolites)
# Purpose:
#   1) To visualize departures from normality across metabolite features,
#      with each point representing one metabolite.
#   2) To identify metabolites that violate either the skewness or kurtosis
#      threshold (highlighted in red).
# Output:
#   outputs/figures/figure1.pdf
# ============================================================

# ---- Libraries ----
library(dplyr)     
library(ggplot2)  
library(moments)  

# ---- File paths (edit to match your repo structure) ----
data_path <- "data/NonPrevExtMaceSamples_IDedDedupMtbs.csv"
ext_mace_ids_path <- "data/ExtMacePatids.txt"
no_ext_mace_ids_path <- "data/NoExtMacePatids.txt"

out_dir <- "outputs/figures"
out_file <- file.path(out_dir, "figure1.pdf")

# ---- Read data ----
dat_all <- read.csv(data_path)

# Patient IDs defining cohorts
ext_mace_patids <- readLines(ext_mace_ids_path)
no_ext_mace_patids <- readLines(no_ext_mace_ids_path)

# ---- Subset (use whichever group you want to assess) ----
# Current choice: no-ExtMACE group
no_ext_mace_subset <- dat_all %>%
  filter(patid %in% no_ext_mace_patids)

# Keep only metabolite columns (drop patid)
no_ext_mace_metabolites <- no_ext_mace_subset %>%
  select(-patid)

# Optional: rename columns to 1:614 (not required for analysis/plotting)
# colnames(no_ext_mace_metabolites) <- as.character(seq_len(ncol(no_ext_mace_metabolites)))

# ---- Compute distribution diagnostics per metabolite ----
skewness_values <- apply(no_ext_mace_metabolites, 2, moments::skewness, na.rm = TRUE)
kurtosis_values <- apply(no_ext_mace_metabolites, 2, moments::kurtosis, na.rm = TRUE)

metrics <- data.frame(
  skewness = skewness_values,
  kurtosis = kurtosis_values
)

# ---- Threshold rule (as in your original code) ----
flag <- (abs(metrics$skewness) > 0.47) & (metrics$kurtosis > 4)

# ---- Plot ----
p <- ggplot(metrics, aes(x = skewness, y = kurtosis)) +
  geom_point(aes(color = flag), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c(`FALSE` = "grey60", `TRUE` = "red"),
    labels = c(`FALSE` = "Within threshold", `TRUE` = "Exceeds thresholds"),
    name = "Normality check"
  ) +
  labs(
    title = "Skewness vs. Kurtosis of Metabolites",
    subtitle = "Red: |skewness| > 0.47 and kurtosis > 4 (computed in no-ExtMACE group)",
    x = "Skewness",
    y = "Kurtosis"
  ) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkred") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ---- Save ----
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(filename = out_file, plot = p, width = 7, height = 7)

print(p)