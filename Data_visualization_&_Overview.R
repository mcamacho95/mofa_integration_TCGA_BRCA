#!/usr/bin/env Rscript
# ==============================================================================
# Script: Data_visualization_&_Overview.R
# Description: Merges clinical annotations, computes UMAP projections for each
#              omics layer, and provides a baseline demographic summary.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(uwot)
})

# Define directories
data_dir <- file.path("data")
fig_dir  <- file.path("results", "figures")
tab_dir  <- file.path("results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/4] Loading Data & Constructing Clinical Master...")
# ------------------------------------------------------------------------------
mofa_inputs <- readRDS(file.path(data_dir, "preprocessed_mofa_inputs.rds"))
clin_tcga   <- readRDS(file.path(data_dir, "clinical_data_tnbc.rds"))
tnbc_class  <- readRDS(file.path(data_dir, "tnbc_classification.rds"))

# Note: Assuming 'samplesheet' is stored in a variable or path (e.g., extracted 
# from the original TCGA colData in Step 01).
# samplesheet <- readRDS(file.path(data_dir, "samplesheet.rds"))

# All samples aligned in step 01 (15-character barcodes)
cohort_samples <- colnames(mofa_inputs$mRNA)

# Build Master Clinical Table using raw variables directly
clinical_master <- tibble(sample = cohort_samples) %>%
  mutate(bcr_patient_barcode = substr(sample, 1, 12)) %>%
  left_join(tnbc_class, by = "bcr_patient_barcode") %>%
  left_join(clin_tcga, by = "bcr_patient_barcode") %>%
  as.data.frame()

saveRDS(clinical_master, file.path(data_dir, "clinical_master.rds"))
message("  -> clinical_master.rds saved successfully.")

message("[2/4] Preparing Matrices for UMAP (Handling NAs)...")
# ------------------------------------------------------------------------------
# Transpose matrices so rows = samples, cols = features
mat_list <- lapply(mofa_inputs, t)

# Drop missing samples entirely for mRNA, lncRNA, miRNA, CNV, Mutations
clean_mats <- list(
  mRNA      = na.omit(mat_list$mRNA),
  lncRNA    = na.omit(mat_list$lncRNA),
  miRNA     = na.omit(mat_list$miRNA),
  CNV       = na.omit(mat_list$CNV),
  Mutations = na.omit(mat_list$Mutations)
)

for(nm in names(clean_mats)) {
  cat(sprintf("  -> %s: Kept %d complete samples across %d features\n", 
              nm, nrow(clean_mats[[nm]]), ncol(clean_mats[[nm]])))
}

# Methylation: Filter features with >50% NAs, then median-impute remainder
meth_mat <- mat
