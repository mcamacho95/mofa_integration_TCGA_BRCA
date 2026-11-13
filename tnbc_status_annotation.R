#!/usr/bin/env Rscript
# ==============================================================================
# Script: tnbc_status_annotation.R
# Description: Retrieves clinical data from TCGA-BRCA, resolves HER2 status
#              using IHC and FISH, and classifies samples as TNBC, Non-TNBC, 
#              or Unclassified.
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(TCGAbiolinks)
})

# Define directories relative to project root
data_dir <- file.path("data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/3] Loading or Fetching TCGA-BRCA Clinical Data...")
# ------------------------------------------------------------------------------
clin_rds_path <- file.path(data_dir, "clinical_data_tnbc.rds")

if (!file.exists(clin_rds_path)) {
  message("  -> Clinical data not found locally. Querying GDC...")
  query_clin <- GDCquery(
    project = "TCGA-BRCA", 
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR XML"
  )
  GDCdownload(query_clin, method = "api", files.per.chunk = 10)
  clinical_data_tnbc <- GDCprepare_clinic(query_clin, clinical.info = "patient")
  saveRDS(clinical_data_tnbc, file = clin_rds_path)
} else {
  message("  -> Loading local clinical data...")
  clinical_data_tnbc <- readRDS(clin_rds_path)
}

message("[2/3] Resolving Receptors and Classifying TNBC Status...")
# ------------------------------------------------------------------------------
# We extract the relevant columns and use dplyr for robust named manipulation
tnbc_df <- clinical_data_tnbc %>%
  dplyr::select(
    bcr_patient_barcode,
    ESR1 = breast_carcinoma_estrogen_receptor_status,
    PGR = breast_carcinoma_progesterone_receptor_status,
    HER2_level = her2_immunohistochemistry_level_result,
    HER2_status = lab_proc_her2_neu_immunohistochemistry_receptor_status,
    HER2_FISH = lab_procedure_her2_neu_in_situ_hybrid_outcome_type
  ) %>%
  # Clean NA and empty strings to "Indeterminate"
  mutate(across(c(ESR1, PGR, HER2_level, HER2_status, HER2_FISH), 
                ~ ifelse(is.na(.) | . == "" | . == "Equivocal", "Indeterminate", .)))

# Resolve HER2 Status
tnbc_df <- tnbc_df %>%
  mutate(
    # 1. Derive HER2 status from IHC level if necessary
    derived_HER2_level = case_when(
      HER2_level %in% c("0", "1+") ~ "Negative",
      HER2_level == "2+" ~ "Indeterminate",
      HER2_level == "3+" ~ "Positive",
      TRUE ~ "Indeterminate"
    ),
    # 2. Prioritize FISH -> Then IHC Status -> Then Derived IHC Level
    final_HER2 = case_when(
      HER2_FISH %in% c("Positive", "Negative") ~ HER2_FISH,
      HER2_status %in% c("Positive", "Negative") ~ HER2_status,
      derived_HER2_level %in% c("Positive", "Negative") ~ derived_HER2_level,
      TRUE ~ "Indeterminate"
    )
  )

# Classify TNBC
tnbc_df <- tnbc_df %>%
  mutate(
    TNBC_status = case_when(
      # If any receptor is explicitly Positive, it is NOT TNBC
      ESR1 == "Positive" | PGR == "Positive" | final_HER2 == "Positive" ~ "Non-TNBC",
      
      # If all three are definitively Negative, it IS TNBC
      ESR1 == "Negative" & PGR == "Negative" & final_HER2 == "Negative" ~ "TNBC",
      
      # Anything else (contains Indeterminates without a Positive) is Unclassified
      TRUE ~ "Unclassified"
    )
  )

# Provide a summary of the classification
message("Classification Summary:")
print(table(tnbc_df$TNBC_status, useNA = "always"))

message("[3/3] Saving Final Annotation...")
# ------------------------------------------------------------------------------
# Keep only the barcode and the final status to merge easily later
final_annotation <- tnbc_df %>%
  dplyr::select(bcr_patient_barcode, TNBC_status)

output_file <- file.path(data_dir, "tnbc_classification.rds")
saveRDS(final_annotation, file = output_file)
message(sprintf("Done! Output saved to: %s", output_file))
