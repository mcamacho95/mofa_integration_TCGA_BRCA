#!/usr/bin/env Rscript
# ==============================================================================
# Script: Data_preprocessing.R
# Description: Pre-processes and normalizes 6 multi-omic modalities (mRNA,
#              lncRNA, miRNA, DNA Methylation, CNV, Mutations) for TCGA-BRCA,
#              harmonizes sample barcodes, and builds aligned input matrices for MOFA+.
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(edgeR)
  library(DESeq2)
  library(SummarizedExperiment)
  library(matrixStats)
})

# Set seed for reproducible transformations
set.seed(142)

# Define directories relative to project root
data_dir <- file.path("data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/7] Processing mRNA and lncRNA (STAR - Counts)...")
# ------------------------------------------------------------------------------
# Expected file: data/data_mrna.rds (SummarizedExperiment from TCGAbiolinks)
mrna_rds_path <- file.path(data_dir, "data_mrna.rds")
if (!file.exists(mrna_rds_path)) {
  stop("Missing 'data_mrna.rds' in data/. Please retrieve via TCGAbiolinks.")
}
data_mrna <- readRDS(mrna_rds_path)

# Extract raw counts and samplesheet
gexp <- as.data.frame(assay(data_mrna, "unstranded")) %>% 
  rownames_to_column("Gene")
samplesheet <- as.data.frame(colData(data_mrna))

# Filter out pseudoautosomal Y genes and problematic indices
genes_table <- as.data.frame(rowRanges(data_mrna)) %>%
  dplyr::slice(-c(15280, 50059, 56397, 57740, 57842, 57950)) %>%
  dplyr::filter(!grepl("_PAR_Y$", gene_id))

# Filter for Primary Solid Tumors (sample type '01')
samplesheet_filtered <- samplesheet %>%
  filter(tissue_type == "Tumor", sample_type == "Primary Tumor")

map_table <- genes_table %>% 
  filter(gene_type %in% c("protein_coding", "lncRNA")) %>%
  dplyr::select(gene_id, gene_name, gene_type)

# Aggregate duplicate gene symbols by summing counts
gexp_aggregated <- gexp %>%
  inner_join(map_table, by = c("Gene" = "gene_id")) %>%
  dplyr::select(-Gene) %>% 
  group_by(gene_name, gene_type) %>% 
  summarise(across(everything(), sum), .groups = 'drop')

colnames(gexp_aggregated) <- gsub("\\.", "-", colnames(gexp_aggregated))
target_samples <- intersect(colnames(gexp_aggregated), rownames(samplesheet_filtered))
gexp_final <- gexp_aggregated %>% dplyr::select(gene_name, gene_type, all_of(target_samples))

feature_types <- gexp_final$gene_type
final_names   <- make.unique(gexp_final$gene_name)
gexp_brca     <- as.matrix(gexp_final %>% dplyr::select(-gene_name, -gene_type))
rownames(gexp_brca) <- final_names

# Standardize to 15-character TCGA barcode (e.g. TCGA-XX-XXXX-01)
colnames(gexp_brca) <- substr(colnames(gexp_brca), 1, 15)

# Deduplicate patients with multiple aliquots by retaining the highest sequencing depth
col_depths <- colSums(gexp_brca, na.rm = TRUE)
gexp_brca  <- gexp_brca[, order(col_depths, decreasing = TRUE)]
gexp_brca  <- gexp_brca[, !duplicated(colnames(gexp_brca))]

# Filter by CPM: protein-coding >= 1 CPM in 20% samples; lncRNA >= 1 CPM in 5% samples
cpms        <- cpm(gexp_brca)
keep_mrna   <- (feature_types == "protein_coding") & (rowSums(cpms >= 1) >= 0.20 * ncol(gexp_brca))
keep_lncrna <- (feature_types == "lncRNA") & (rowSums(cpms >= 1) >= 0.05 * ncol(gexp_brca))
keep        <- keep_mrna | keep_lncrna

gexp_filtered          <- gexp_brca[keep, ]
feature_types_filtered <- feature_types[keep]

# Variance-stabilizing transformation (DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = gexp_filtered,
  colData   = data.frame(row.names = colnames(gexp_filtered), sample = colnames(gexp_filtered)),
  design    = ~ 1
)
dds      <- estimateSizeFactors(dds)
vst_out  <- vst(dds, blind = TRUE)
gexp_vst <- assay(vst_out)

mRNA_vst   <- gexp_vst[feature_types_filtered == "protein_coding", ]
lncRNA_vst <- gexp_vst[feature_types_filtered == "lncRNA", ]
message(sprintf("  -> Retained %d mRNA and %d lncRNA features across %d samples.", 
                nrow(mRNA_vst), nrow(lncRNA_vst), ncol(mRNA_vst)))


message("[2/7] Processing miRNA (BCGSC Profiling)...")
# ------------------------------------------------------------------------------
mirna_rds_path <- file.path(data_dir, "data_mirna_brca.rds")
if (!file.exists(mirna_rds_path)) {
  stop("Missing 'data_mirna_brca.rds' in data/.")
}
data_mirna <- readRDS(mirna_rds_path)

miRNA_brca <- data_mirna %>%
  as.data.frame() %>%
  column_to_rownames("miRNA_ID") %>%
  dplyr::select(starts_with("read_count"))

colnames(miRNA_brca) <- gsub("read_count_", "", colnames(miRNA_brca))
colnames(miRNA_brca) <- gsub("\\.", "-", colnames(miRNA_brca))

# Retain primary tumors and standardize to 15-char barcode
primary_mir_cols <- colnames(miRNA_brca)[substr(colnames(miRNA_brca), 14, 15) == "01"]
miRNA_brca <- miRNA_brca[, primary_mir_cols]
colnames(miRNA_brca) <- substr(colnames(miRNA_brca), 1, 15)

# Deduplicate by highest total counts
sample_depths_mir <- colSums(miRNA_brca, na.rm = TRUE)
miRNA_brca        <- miRNA_brca[, order(sample_depths_mir, decreasing = TRUE)]
miRNA_brca        <- miRNA_brca[, !duplicated(colnames(miRNA_brca))]

# Filter: >= 1 CPM in 5% of samples
cpms_mir   <- cpm(miRNA_brca)
keep_mir   <- rowSums(cpms_mir >= 1) >= 0.05 * ncol(miRNA_brca)
miRNA_filt <- miRNA_brca[keep_mir, ]

# VST Normalization
dds_mir   <- DESeqDataSetFromMatrix(
  countData = miRNA_filt,
  colData   = data.frame(row.names = colnames(miRNA_filt), sample = colnames(miRNA_filt)),
  design    = ~ 1
)
dds_mir   <- estimateSizeFactors(dds_mir)
vsd_mir   <- varianceStabilizingTransformation(dds_mir, blind = TRUE)
miRNA_vst <- assay(vsd_mir)
message(sprintf("  -> Retained %d miRNA features across %d samples.", nrow(miRNA_vst), ncol(miRNA_vst)))


message("[3/7] Processing Somatic Mutations...")
# ------------------------------------------------------------------------------
mut_file <- file.path(data_dir, "mut_brca.csv")
mutations <- read.table(mut_file, sep = ",", header = TRUE, quote = "", 
                        row.names = 1, check.names = FALSE, fill = TRUE)
colnames(mutations) <- gsub("\\.", "-", colnames(mutations))
mut_tumor_cols <- colnames(mutations)[substr(colnames(mutations), 14, 15) == "01"]
mutations <- mutations[, mut_tumor_cols]
colnames(mutations) <- substr(colnames(mutations), 1, 15)
mutations <- mutations[, !duplicated(colnames(mutations))]

# Retain genes mutated in at least 10 samples
mutations <- mutations[rowSums(mutations > 0) >= 10, ]
mutations <- as.matrix(mutations)
message(sprintf("  -> Retained %d mutated genes across %d samples.", nrow(mutations), ncol(mutations)))


message("[4/7] Processing Copy Number Variation (CNV)...")
# ------------------------------------------------------------------------------
cnv_file <- file.path(data_dir, "cnv_brca.csv")
CNV <- read.table(cnv_file, sep = ",", header = TRUE, quote = "", 
                  row.names = 1, check.names = FALSE, fill = TRUE)
colnames(CNV) <- gsub("\\.", "-", colnames(CNV))
cnv_tumor_cols <- colnames(CNV)[substr(colnames(CNV), 14, 15) == "01"]
CNV <- CNV[, cnv_tumor_cols]
colnames(CNV) <- substr(colnames(CNV), 1, 15)
CNV <- CNV[, !duplicated(colnames(CNV))]

# Standardize features across samples (Z-score scaling per gene)
CNV_scaled <- t(scale(t(as.matrix(CNV))))
is_nan     <- apply(CNV_scaled, 1, function(x) any(is.na(x)))
CNV_scaled <- CNV_scaled[!is_nan, ]
message(sprintf("  -> Retained %d CNV features across %d samples.", nrow(CNV_scaled), ncol(CNV_scaled)))


message("[5/7] Processing DNA Methylation (M-values)...")
# ------------------------------------------------------------------------------
meth_file <- file.path(data_dir, "meth.csv")
meth <- read.table(meth_file, sep = ",", header = TRUE, quote = "", 
                   row.names = 1, check.names = FALSE, fill = TRUE)
colnames(meth) <- gsub('"', '', colnames(meth))
rownames(meth) <- gsub('"', '', rownames(meth))
colnames(meth) <- gsub("\\.", "-", colnames(meth))

meth_tumor_cols <- colnames(meth)[substr(colnames(meth), 14, 15) == "01"]
meth <- meth[, meth_tumor_cols]
colnames(meth) <- substr(colnames(meth), 1, 15)
meth <- meth[, !duplicated(colnames(meth))]

# Convert Beta-values to M-values with offset clipping
offset <- 0.001
meth[meth < offset]      <- offset
meth[meth > 1 - offset]  <- 1 - offset
meth_mvalues <- log2(as.matrix(meth) / (1 - as.matrix(meth)))
message(sprintf("  -> Retained %d methylation features across %d samples.", 
                nrow(meth_mvalues), ncol(meth_mvalues)))


message("[6/7] Multi-Omic Alignment across cohort...")
# ------------------------------------------------------------------------------
# Union of all unique primary tumor samples
all_samples <- sort(unique(c(
  colnames(mRNA_vst), 
  colnames(lncRNA_vst), 
  colnames(miRNA_vst),
  colnames(mutations), 
  colnames(CNV_scaled), 
  colnames(meth_mvalues)
)))

# Helper function to align matrix and fill missing samples with NA
complete_matrix <- function(mat, cohort_samples) {
  missing_samples <- setdiff(cohort_samples, colnames(mat))
  if (length(missing_samples) > 0) {
    na_mat <- matrix(NA, nrow = nrow(mat), ncol = length(missing_samples))
    colnames(na_mat) <- missing_samples
    rownames(na_mat) <- rownames(mat)
    full_mat <- cbind(mat, na_mat)
  } else {
    full_mat <- mat
  }
  return(full_mat[, cohort_samples, drop = FALSE])
}

mofa_input_list <- list(
  mRNA        = complete_matrix(mRNA_vst, all_samples),
  lncRNA      = complete_matrix(lncRNA_vst, all_samples),
  miRNA       = complete_matrix(miRNA_vst, all_samples),
  Methylation = complete_matrix(meth_mvalues, all_samples),
  CNV         = complete_matrix(CNV_scaled, all_samples),
  Mutations   = complete_matrix(mutations, all_samples)
)

message(sprintf("Successfully aligned all 6 views across %d total unique patients.", length(all_samples)))


message("[7/7] Saving Aligned Data Object...")
# ------------------------------------------------------------------------------
output_file <- file.path(data_dir, "preprocessed_mofa_inputs.rds")
saveRDS(mofa_input_list, file = output_file)
message(sprintf("Done! Output saved to: %s", output_file))
