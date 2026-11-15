#!/usr/bin/env Rscript
# ==============================================================================
# Script: GSEA.R
# Description: Performs Gene Set Enrichment Analysis (fgsea) on MOFA+ feature
#              loadings against MSigDB Hallmark and GO:BP collections for:
#                1) Factor 1 mRNA loadings (two-sided enrichment: Luminal vs. TNBC)
#                2) Factor 2 Methylation loadings (one-sided positive enrichment)
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(tibble)
})

# Define directories
model_dir <- file.path("results", "models")
tab_dir   <- file.path("results", "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/4] Loading Optimal MOFA+ Model...")
# ------------------------------------------------------------------------------
model_path <- file.path(model_dir, "best_mofa_model.hdf5")
if (!file.exists(model_path)) {
  stop("Missing 'best_mofa_model.hdf5' in results/models/. Run 04_mofa_train.R first.")
}
model_final <- load_model(model_path)

message("[2/4] Retrieving MSigDB Gene Sets (Hallmark & GO:BP)...")
# ------------------------------------------------------------------------------
# Hallmark gene sets (Category H)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  distinct(gs_name, gene_symbol)
pathways_hallmark <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Gene Ontology: Biological Process (Category C5, Subcategory GO:BP)
gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  distinct(gs_name, gene_symbol)
pathways_gobp <- split(gobp_sets$gene_symbol, gobp_sets$gs_name)

message(sprintf("  -> Loaded %d Hallmark and %d GO:BP gene sets.",
                length(pathways_hallmark), length(pathways_gobp)))

# Helper Functions
# ------------------------------------------------------------------------------
clean_pathway_names <- function(x) {
  x <- gsub("^(HALLMARK_|GOBP_)", "", x)
  x <- gsub("_", " ", x)
  return(x)
}

format_fdr <- function(fdr_vec) {
  case_when(
    fdr_vec < 0.001 ~ "< 0.001",
    TRUE            ~ sprintf("%.3f", fdr_vec)
  )
}

run_factor_gsea <- function(model, view_name, factor_idx, pathways, collection_label,
                            score_type = "std", seed = 8, min_size = 15, max_size = 500) {
  message(sprintf("\n--- Running GSEA: Factor %d | View: %s | Collection: %s (scoreType = '%s') ---",
                  factor_idx, view_name, collection_label, score_type))
  
  # 1. Extract and rank feature weights in descending order
  loadings_mat <- get_weights(model, views = view_name, factors = factor_idx)[[1]]
  ranks        <- loadings_mat[, 1]
  names(ranks) <- gsub(paste0("_", view_name, "$"), "", rownames(loadings_mat))
  ranks        <- ranks[!is.na(ranks)]
  ranks        <- sort(ranks, decreasing = TRUE)
  
  # 2. Run fgsea with single-threaded execution (nproc = 1) for exact seed reproducibility
  set.seed(seed)
  fgsea_res <- fgsea(
    pathways  = pathways,
    stats     = ranks,
    scoreType = score_type,
    minSize   = min_size,
    maxSize   = max_size,
    nproc     = 1
  ) %>%
    as_tibble() %>%
    filter(padj < 0.05) %>%
    arrange(desc(NES))
  
  message(sprintf("  -> Significant pathways (FDR < 0.05): %d", nrow(fgsea_res)))
  
  if (nrow(fgsea_res) == 0) {
    return(invisible(NULL))
  }
  
  # 3. Export full results (flattening leadingEdge list column for CSV compatibility)
  fgsea_export <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
  
  full_file <- file.path(tab_dir, sprintf("gsea_Factor%d_%s_%s_full.csv",
                                          factor_idx, view_name, collection_label))
  write.csv(fgsea_export, full_file, row.names = FALSE)
  
  # 4. Format Top 10 Positively and Negatively Enriched Pathways for Thesis Tables
  top_pos <- fgsea_res %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    slice_head(n = 10) %>%
    mutate(
      Direction = "Positive",
      Pathway   = clean_pathway_names(pathway),
      NES       = round(NES, 2),
      FDR       = format_fdr(padj)
    ) %>%
    dplyr::select(Direction, Pathway, NES, FDR, Size = size)
  
  if (score_type == "std") {
    top_neg <- fgsea_res %>%
      filter(NES < 0) %>%
      arrange(NES) %>%
      slice_head(n = 10) %>%
      mutate(
        Direction = "Negative",
        Pathway   = clean_pathway_names(pathway),
        NES       = round(NES, 2),
        FDR       = format_fdr(padj)
      ) %>%
      dplyr::select(Direction, Pathway, NES, FDR, Size = size)
    
    formatted_summary <- bind_rows(top_pos, top_neg)
  } else {
    formatted_summary <- top_pos
  }
  
  summary_file <- file.path(tab_dir, sprintf("gsea_Factor%d_%s_%s_top10_formatted.csv",
                                             factor_idx, view_name, collection_label))
  write.csv(formatted_summary, summary_file, row.names = FALSE)
  print(as.data.frame(formatted_summary))
  
  return(invisible(fgsea_res))
}

message("[3/4] Running Two-Sided GSEA on Factor 1 (mRNA View)...")
# ------------------------------------------------------------------------------
run_factor_gsea(
  model            = model_final,
  view_name        = "mRNA",
  factor_idx       = 1,
  pathways         = pathways_hallmark,
  collection_label = "Hallmark",
  score_type       = "std",
  seed             = 8
)

run_factor_gsea(
  model            = model_final,
  view_name        = "mRNA",
  factor_idx       = 1,
  pathways         = pathways_gobp,
  collection_label = "GOBP",
  score_type       = "std",
  seed             = 9
)

message("[4/4] Running One-Sided Positive GSEA on Factor 2 (Methylation View)...")
# ------------------------------------------------------------------------------
# Over 99% of Factor 2 Methylation weights are positive; assess positive pole only
run_factor_gsea(
  model            = model_final,
  view_name        = "Methylation",
  factor_idx       = 2,
  pathways         = pathways_hallmark,
  collection_label = "Hallmark",
  score_type       = "pos",
  seed             = 8
)

run_factor_gsea(
  model            = model_final,
  view_name        = "Methylation",
  factor_idx       = 2,
  pathways         = pathways_gobp,
  collection_label = "GOBP",
  score_type       = "pos",
  seed             = 9
)

message("\nDone! All GSEA tables saved to results/tables/.")
