#!/usr/bin/env Rscript
# ==============================================================================
# Script: DGE_limma.R
# Description: Performs supervised differential expression analysis (limma)
#              comparing TNBC vs. Non-TNBC tumours for mRNA and lncRNA views,
#              validates concordance with MOFA+ Factor 1 weights using an
#              unbroken FDR < 0.05 rank cutoff, and exports volcano plots and
#              prioritized TNBC gene lists for downstream drug repurposing.
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(limma)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# Define directories
data_dir  <- file.path("data")
model_dir <- file.path("results", "models")
fig_dir   <- file.path("results", "figures")
tab_dir   <- file.path("results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/3] Loading Optimal MOFA+ Model and Clinical Metadata...")
# ------------------------------------------------------------------------------
model_path <- file.path(model_dir, "best_mofa_model.hdf5")
clin_path  <- file.path(data_dir, "clinical_master.rds")

if (!file.exists(model_path)) stop("Missing 'best_mofa_model.hdf5'. Run 04_mofa_train.R first.")
if (!file.exists(clin_path))  stop("Missing 'clinical_master.rds'. Run 03_data_visualization.R first.")

model_final     <- load_model(model_path)
clinical_master <- readRDS(clin_path)

# Align clinical metadata to MOFA+ sample ordering and filter to confirmed TNBC/Non-TNBC
model_samples <- samples_names(model_final)[[1]]
meta_df <- clinical_master[match(model_samples, clinical_master$sample), ] %>%
  mutate(
    tnbc = case_when(
      TNBC_status %in% c("TNBC") ~ "TNBC",
      TNBC_status %in% c("Non-TNBC", "NO TNBC") ~ "Non-TNBC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(tnbc))

message(sprintf("  -> Retained %d samples with confirmed TNBC status (%d TNBC, %d Non-TNBC).",
                nrow(meta_df),
                sum(meta_df$tnbc == "TNBC"),
                sum(meta_df$tnbc == "Non-TNBC")))

# Modular DGE + Factor 1 Concordance + Volcano Plot Function
# ------------------------------------------------------------------------------
run_view_dge <- function(view_name, model, metadata, plot_title, plot_filename,
                         logFC_cut = 1.0, fdr_cut = 0.05) {
  message(sprintf("\n--- Running limma DGE & Factor 1 Validation for %s ---", view_name))
  
  # 1. Extract expression matrix and strip any view suffix
  expr_matrix <- get_data(model, views = view_name)[[1]][[1]]
  rownames(expr_matrix) <- sub(paste0("_", view_name, "$"), "", rownames(expr_matrix))
  
  # 2. Synchronize matrix columns to metadata samples
  valid_samples     <- intersect(metadata$sample, colnames(expr_matrix))
  meta_sync         <- metadata %>% filter(sample %in% valid_samples)
  expr_matrix_clean <- expr_matrix[, meta_sync$sample, drop = FALSE]
  
  # 3. Fit limma linear model with empirical Bayes moderation
  group  <- factor(meta_sync$tnbc, levels = c("Non-TNBC", "TNBC"))
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", "TNBC_vs_NonTNBC")
  
  fit         <- lmFit(expr_matrix_clean, design)
  fit         <- eBayes(fit)
  dea_results <- topTable(fit, coef = "TNBC_vs_NonTNBC", number = Inf, sort.by = "P")
  
  write.csv(
    dea_results %>% rownames_to_column("Feature"),
    file.path(tab_dir, sprintf("dge_limma_%s_full.csv", view_name)),
    row.names = FALSE
  )
  
  # 4. Merge with Factor 1 weights and apply unbroken FDR < 0.05 rank cutoff
  f1_weights <- get_weights(model, views = view_name, factors = 1, as.data.frame = TRUE) %>%
    mutate(feature = sub(paste0("_", view_name, "$"), "", feature)) %>%
    dplyr::select(feature, MOFA_weight = value)
  
  validation_df <- dea_results %>%
    rownames_to_column("feature") %>%
    inner_join(f1_weights, by = "feature") %>%
    mutate(
      abs_weight  = abs(MOFA_weight),
      significant = adj.P.Val < fdr_cut
    ) %>%
    arrange(desc(abs_weight)) %>%
    mutate(rank = row_number())
  
  first_insig_rank <- validation_df %>%
    filter(!significant) %>%
    slice_head(n = 1) %>%
    pull(rank)
  
  # Handle edge case where all features are significant
  if (length(first_insig_rank) == 0) first_insig_rank <- nrow(validation_df) + 1
  
  # Negative Factor 1 weights correspond to the Basal/TNBC pole
  validation_subset <- validation_df %>%
    filter(rank < first_insig_rank) %>%
    mutate(
      MOFA_pole = if_else(MOFA_weight < 0, "TNBC_enriched", "NonTNBC_enriched"),
      DGE_direction = if_else(logFC > 0, "TNBC_up", "NonTNBC_up"),
      concordant = if_else(
        (MOFA_pole == "TNBC_enriched" & DGE_direction == "TNBC_up") |
        (MOFA_pole == "NonTNBC_enriched" & DGE_direction == "NonTNBC_up"),
        "Concordant", "Discordant"
      )
    )
  
  cat(sprintf("  Unbroken Window (Ranks 1 to %d): %d features\n",
              first_insig_rank - 1, nrow(validation_subset)))
  cat(sprintf("  Concordant: %d (%.2f%%) | Discordant: %d\n",
              sum(validation_subset$concordant == "Concordant"),
              mean(validation_subset$concordant == "Concordant") * 100,
              sum(validation_subset$concordant == "Discordant")))
  print(table(MOFA_Pole = validation_subset$MOFA_pole,
              DGE_Direction = validation_subset$DGE_direction))
  
  write.csv(
    validation_subset,
    file.path(tab_dir, sprintf("dge_factor1_concordant_%s.csv", view_name)),
    row.names = FALSE
  )
  
  # Save TNBC-upregulated concordant features (used for DGIdb in Step 09)
  tnbc_up_features <- validation_subset %>%
    filter(concordant == "Concordant", DGE_direction == "TNBC_up") %>%
    pull(feature)
  
  writeLines(
    paste(tnbc_up_features, collapse = ", "),
    file.path(tab_dir, sprintf("tnbc_upregulated_%s_dgidb_input.txt", view_name))
  )
  saveRDS(
    tnbc_up_features,
    file.path(tab_dir, sprintf("tnbc_upregulated_%s.rds", view_name))
  )
  
  # 5. Generate and Save Volcano Plot
  volcano_df <- dea_results %>%
    rownames_to_column(var = "Gene") %>%
    mutate(
      neg_log10_adjP = -log10(adj.P.Val),
      Significance = case_when(
        adj.P.Val < fdr_cut & logFC >  logFC_cut ~ "Upregulated in TNBC",
        adj.P.Val < fdr_cut & logFC < -logFC_cut ~ "Downregulated in TNBC",
        TRUE                                    ~ "Not Significant"
      )
    )
  
  top_up <- volcano_df %>%
    filter(Significance == "Upregulated in TNBC") %>%
    slice_max(order_by = abs(logFC) * neg_log10_adjP, n = 10)
  
  top_down <- volcano_df %>%
    filter(Significance == "Downregulated in TNBC") %>%
    slice_max(order_by = abs(logFC) * neg_log10_adjP, n = 10)
  
  top_genes <- bind_rows(top_up, top_down)
  
  p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_adjP)) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.6) +
    scale_color_manual(values = c(
      "Upregulated in TNBC"   = "#D95F02",
      "Downregulated in TNBC" = "#1F78B4",
      "Not Significant"       = "#BDBDBD"
    )) +
    geom_vline(xintercept = c(-logFC_cut, logFC_cut), linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene),
      size = 3.2,
      box.padding = 0.35,
      point.padding = 0.3,
      max.overlaps = 20,
      segment.color = "grey50",
      seed = 42
    ) +
    labs(
      title    = plot_title,
      subtitle = paste0("FDR < ", fdr_cut, ", |log2FC| > ", logFC_cut),
      x        = expression(Log[2]~"Fold Change (TNBC / Non-TNBC)"),
      y        = expression(-Log[10]~"Adjusted "*italic(p)*"-value"),
      color    = "Regulation"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(fig_dir, plot_filename), plot = p_volcano, width = 7, height = 5.5, dpi = 300)
  message(sprintf("  -> Saved volcano plot to results/figures/%s", plot_filename))
  
  return(invisible(validation_subset))
}

message("[2/3] Running DGE Analysis on mRNA View...")
# ------------------------------------------------------------------------------
run_view_dge(
  view_name     = "mRNA",
  model         = model_final,
  metadata      = meta_df,
  plot_title    = "Differential mRNA Expression: TNBC vs Non-TNBC",
  plot_filename = "DGE_Factor1_mRNAs.png"
)

message("[3/3] Running DGE Analysis on lncRNA View...")
# ------------------------------------------------------------------------------
run_view_dge(
  view_name     = "lncRNA",
  model         = model_final,
  metadata      = meta_df,
  plot_title    = "Differential lncRNA Expression: TNBC vs Non-TNBC",
  plot_filename = "DGE_Factor1_lncRNA.png"
)

message("\nDone! All DGE results and gene lists exported to results/.")
