#!/usr/bin/env Rscript
# ==============================================================================
# Script: Drug_repurposing.R
# Description: Integrates Factor 1 TNBC-upregulated mRNA biomarkers (from Step 06)
#              with Drug-Gene Interaction Database (DGIdb) annotations to:
#                1) Quantify approved vs. investigational drug-gene interactions
#                2) Identify highly connected druggable target hubs (> Q75)
#                3) Prioritize targets combining Factor 1 loadings & DGIdb scores
#                4) Identify multi-target approved compounds across TNBC drivers
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# Define directories
data_dir <- file.path("data")
tab_dir  <- file.path("results", "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/5] Loading Factor 1 TNBC-Upregulated mRNA Candidates & DGIdb Export...")
# ------------------------------------------------------------------------------
dge_path   <- file.path(tab_dir, "dge_factor1_concordant_mRNA.csv")
dgidb_path <- file.path(data_dir, "gene_interactions.tsv")

if (!file.exists(dge_path)) {
  stop("Missing 'dge_factor1_concordant_mRNA.csv' in results/tables/. Run 06_dge_limma.R first.")
}

if (!file.exists(dgidb_path)) {
  stop(paste0(
    "Missing 'gene_interactions.tsv' in data/.\n",
    "To generate this file:\n",
    "  1. Copy the gene list from 'results/tables/tnbc_upregulated_mRNA_dgidb_input.txt'.\n",
    "  2. Query the Drug-Gene Interaction Database (https://www.dgidb.org/).\n",
    "  3. Download the TSV interactions table and save it as 'data/gene_interactions.tsv'."
  ))
}

validation_subset <- read_csv(dge_path, show_col_types = FALSE)

# Extract the 677 concordant TNBC-upregulated protein-coding genes
val_clean <- validation_subset %>%
  filter(concordant == "Concordant", DGE_direction == "TNBC_up") %>%
  mutate(gene = gsub("_mRNA$", "", feature)) %>%
  dplyr::select(gene, MOFA_weight, abs_weight, logFC, adj.P.Val, MOFA_pole, DGE_direction)

tnbc_genes  <- unique(val_clean$gene)
n_submitted <- length(tnbc_genes)

# Read DGIdb TSV export and standardize column names (spaces -> underscores, lowercase)
dgi_raw <- read_tsv(dgidb_path, show_col_types = FALSE)
dgi_df  <- dgi_raw %>%
  rename_with(~ tolower(gsub(" ", "_", .x)))

# Standardize gene/drug column names across DGIdb export versions
if (!"gene" %in% colnames(dgi_df) && "gene_name" %in% colnames(dgi_df)) {
  dgi_df <- dgi_df %>% rename(gene = gene_name)
}
if (!"drug" %in% colnames(dgi_df) && "drug_name" %in% colnames(dgi_df)) {
  dgi_df <- dgi_df %>% rename(drug = drug_name)
}

message("[2/5] Summarizing Overall & Clinically Approved DGIdb Interactions...")
# ------------------------------------------------------------------------------
g_matched_genes  <- n_distinct(dgi_df$gene, na.rm = TRUE)
d_distinct_drugs <- n_distinct(dgi_df$drug, na.rm = TRUE)

approved_clean <- dgi_df %>%
  filter(!is.na(gene), !is.na(drug), regulatory_approval == "Approved")

x_approved_drugs          <- n_distinct(approved_clean$drug)
genes_with_approved_drugs <- n_distinct(approved_clean$gene)

cat(sprintf(
  "\nSummary: Of the %d TNBC-upregulated candidate genes submitted, %d had registered interactions in DGIdb across %d distinct drugs. Of these, %d drugs are clinically approved (targeting %d distinct candidate genes).\n\n",
  n_submitted, g_matched_genes, d_distinct_drugs, x_approved_drugs, genes_with_approved_drugs
))

message("[3/5] Identifying Highly Connected Druggable Target Hubs (> Q75)...")
# ------------------------------------------------------------------------------
gene_drug_pairs <- dgi_df %>%
  filter(!is.na(gene), !is.na(drug)) %>%
  group_by(gene, drug) %>%
  summarise(
    is_approved = any(regulatory_approval == "Approved", na.rm = TRUE),
    .groups     = "drop"
  )

summary_df <- gene_drug_pairs %>%
  group_by(gene) %>%
  summarise(
    total_interactions        = n_distinct(drug),
    approved_interactions     = sum(is_approved),
    non_approved_interactions = sum(!is_approved),
    .groups                   = "drop"
  ) %>%
  arrange(desc(approved_interactions), desc(total_interactions))

write_csv(summary_df, file.path(tab_dir, "dgidb_gene_interaction_summary.csv"))

# Filter to genes with >= 1 approved drug and compute distribution quantiles
approved_targets <- summary_df %>%
  filter(approved_interactions > 0)

med_int <- median(approved_targets$approved_interactions)
q25_int <- quantile(approved_targets$approved_interactions, 0.25)[[1]]
q75_int <- quantile(approved_targets$approved_interactions, 0.75)[[1]]
max_int <- max(approved_targets$approved_interactions)

message(sprintf("  -> Approved drugs per target gene (n = %d): Median = %.1f (IQR: %.1f - %.1f, Max = %d)",
                nrow(approved_targets), med_int, q25_int, q75_int, max_int))

hub_genes <- approved_targets %>%
  filter(approved_interactions > q75_int) %>%
  arrange(desc(approved_interactions))

write_csv(hub_genes, file.path(tab_dir, "dgidb_approved_hub_genes.csv"))
message(sprintf("  -> Exported %d target hub genes (> Q75) to results/tables/dgidb_approved_hub_genes.csv",
                nrow(hub_genes)))

message("[4/5] Prioritizing Targets by Factor 1 Loadings & DGIdb Interaction Scores...")
# ------------------------------------------------------------------------------
merged_drugs_df <- approved_clean %>%
  inner_join(val_clean, by = "gene")

# Save full ranked drug-gene interaction table for inspection (e.g., filtering lab buffer artifacts)
ranked_interactions <- merged_drugs_df %>%
  arrange(desc(interaction_score), desc(abs_weight))

write_csv(ranked_interactions, file.path(tab_dir, "dgidb_approved_interactions_ranked.csv"))

# Gene-level prioritization combining MOFA+ Factor 1 weight and max DGIdb score
funnel_summary <- merged_drugs_df %>%
  group_by(gene) %>%
  summarise(
    MOFA_weight        = first(MOFA_weight),
    abs_weight         = first(abs_weight),
    logFC              = first(logFC),
    adj_P_Val          = first(adj.P.Val),
    max_dgidb_score    = max(interaction_score, na.rm = TRUE),
    n_approved_drugs   = n_distinct(drug),
    top_approved_drugs = paste(head(unique(drug[order(-interaction_score)]), 5), collapse = "; "),
    bio_evidence_score = first(abs_weight) * max(interaction_score, na.rm = TRUE),
    .groups            = "drop"
  ) %>%
  arrange(desc(bio_evidence_score), desc(n_approved_drugs))

write_csv(funnel_summary, file.path(tab_dir, "dgidb_prioritized_targets_ranked.csv"))

message("[5/5] Identifying Multi-Target Approved Compounds (>= 2 Target Genes)...")
# ------------------------------------------------------------------------------
has_indication <- "indication" %in% colnames(merged_drugs_df)

# A. Multi-target drugs across top 25 prioritized target hubs
top_hub_genes <- head(funnel_summary$gene, 25)

multi_target_top25 <- merged_drugs_df %>%
  filter(gene %in% top_hub_genes) %>%
  group_by(drug) %>%
  summarise(
    n_targets        = n_distinct(gene),
    target_list      = paste(sort(unique(gene)), collapse = ", "),
    max_dgidb_score  = max(interaction_score, na.rm = TRUE),
    primary_indication = if (has_indication) paste(head(unique(na.omit(indication)), 1), collapse = "; ") else NA_character_,
    .groups          = "drop"
  ) %>%
  filter(n_targets >= 2) %>%
  arrange(desc(n_targets), desc(max_dgidb_score))

write_csv(multi_target_top25, file.path(tab_dir, "dgidb_multi_target_drugs_top25_hubs.csv"))

# B. Multi-target drugs across all 150 approved-target genes
multi_target_all <- merged_drugs_df %>%
  group_by(drug) %>%
  summarise(
    n_targets          = n_distinct(gene),
    target_list        = paste(sort(unique(gene)), collapse = ", "),
    max_dgidb_score    = max(interaction_score, na.rm = TRUE),
    mean_dgidb_score   = round(mean(interaction_score, na.rm = TRUE), 3),
    primary_indication = if (has_indication) paste(head(unique(na.omit(indication)), 1), collapse = "; ") else NA_character_,
    .groups            = "drop"
  ) %>%
  filter(n_targets >= 2) %>%
  arrange(desc(n_targets), desc(max_dgidb_score))

write_csv(multi_target_all, file.path(tab_dir, "dgidb_multi_target_drugs_all_approved.csv"))
message(sprintf("  -> Identified %d multi-target approved drugs across all %d targetable TNBC genes.",
                nrow(multi_target_all), genes_with_approved_drugs))

# Optional: Check survival-associated biomarkers if secondary DGIdb export exists
surv_dgidb_path <- file.path(data_dir, "gene_interactions_cox_glmnet.tsv")
if (file.exists(surv_dgidb_path)) {
  message("\n[Optional] Processing Survival Biomarker DGIdb Interactions...")
  dgi_surv <- read_tsv(surv_dgidb_path, show_col_types = FALSE) %>%
    rename_with(~ tolower(gsub(" ", "_", .x)))
  approved_surv <- dgi_surv %>%
    filter(regulatory_approval == "Approved")
  write_csv(approved_surv, file.path(tab_dir, "dgidb_survival_biomarkers_approved.csv"))
  message(sprintf("  -> Exported %d approved interactions for survival biomarkers.", nrow(approved_surv)))
}

message("\nDone! All drug repurposing tables saved to results/tables/.")
