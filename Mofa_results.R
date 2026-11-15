#!/usr/bin/env Rscript
# ==============================================================================
# Script: Mofa_results.R
# Description: Evaluates the optimal trained MOFA+ model: variance explained,
#              clinical phenotype associations, Factor 1 subtype discrimination
#              (violin plots & ROC curves), top feature loadings, and heatmaps.
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(plotROC)
})

# Define directories
data_dir  <- file.path("data")
model_dir <- file.path("results", "models")
fig_dir   <- file.path("results", "figures")
tab_dir   <- file.path("results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/5] Loading Optimal MOFA+ Model and Clinical Metadata...")
# ------------------------------------------------------------------------------
model_path <- file.path(model_dir, "best_mofa_model.hdf5")
clin_path  <- file.path(data_dir, "clinical_master.rds")

if (!file.exists(model_path)) stop("Missing 'best_mofa_model.hdf5'. Run 04_mofa_train.R first.")
if (!file.exists(clin_path))  stop("Missing 'clinical_master.rds'. Run 03_data_visualization.R first.")

model_final     <- load_model(model_path)
clinical_master <- readRDS(clin_path)

# Clean feature suffixes if present (e.g., _mRNA, _CNV)
features_names(model_final) <- lapply(features_names(model_final), function(feats) {
  sub("_(mRNA|CNV|Methylation|lncRNA|miRNA|Mutations)$", "", feats)
})

# Align clinical metadata to exact sample order inside MOFA2
model_samples <- samples_names(model_final)[[1]]
meta_aligned  <- clinical_master[match(model_samples, clinical_master$sample), ]

# Standardize key subtype columns needed for downstream visualizations
meta_aligned <- meta_aligned %>%
  mutate(
    tnbc = case_when(
      TNBC_status %in% c("TNBC") ~ "TNBC",
      TNBC_status %in% c("Non-TNBC", "NO TNBC") ~ "Non-TNBC",
      TRUE ~ NA_character_
    ),
    PAM50 = if ("paper_BRCA_Subtype_PAM50" %in% colnames(.)) {
      as.character(paper_BRCA_Subtype_PAM50)
    } else if ("PAM50" %in% colnames(.)) {
      as.character(PAM50)
    } else {
      NA_character_
    },
    Basal_like = case_when(
      PAM50 == "Basal" ~ "Basal",
      !is.na(PAM50) ~ "Non-Basal",
      TRUE ~ NA_character_
    )
  )

samples_metadata(model_final) <- meta_aligned

message("[2/5] Quantifying Variance Explained...")
# ------------------------------------------------------------------------------
r2_res <- get_variance_explained(model_final)
write.csv(r2_res$r2_total[[1]], file.path(tab_dir, "mofa_total_variance_explained.csv"), row.names = TRUE)
write.csv(r2_res$r2_per_factor[[1]], file.path(tab_dir, "mofa_variance_per_factor.csv"), row.names = TRUE)

p_var_total <- plot_variance_explained(model_final, x = "view", y = "factor", plot_total = TRUE)[[2]]
ggsave(file.path(fig_dir, "Total_Var_explained.png"), plot = p_var_total, width = 6, height = 4.5, dpi = 300)

p_var_factor <- plot_variance_explained(model_final, x = "view", y = "factor")
ggsave(file.path(fig_dir, "Var_explained_per_factor.png"), plot = p_var_factor, width = 6.5, height = 5.5, dpi = 300)
message("  -> Variance explained tables and plots saved.")

message("[3/5] Testing Statistical Associations with Clinical Variables...")
# ------------------------------------------------------------------------------
n_factors   <- get_dimensions(model_final)$K
factors_mat <- get_factors(model_final, factors = 1:n_factors)[[1]]

# Aggregate factor scores at patient level (12-char barcode)
factors_long <- as.data.frame(factors_mat) %>%
  rownames_to_column("sample") %>%
  mutate(patient_key = substr(sample, 1, 12)) %>%
  group_by(patient_key) %>%
  summarise(across(starts_with("Factor"), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = starts_with("Factor"), names_to = "Factor", values_to = "factor_score")

# Exclude technical IDs and time-to-event columns (evaluated separately in survival)
ignore_cols <- c(
  "sample", "patient_key", "barcode", "patient", "sample_submitter_id", "paper_patient",
  "bcr_patient_barcode", "submitter_id", "raw_id", "patient_id", "project",
  "days_to_death", "days_to_last_follow_up", "days_to_death_unified",
  "days_to_followup_unified", "serial_time", "vital_status_unified"
)

clinical_patient <- meta_aligned %>%
  mutate(patient_key = substr(sample, 1, 12)) %>%
  distinct(patient_key, .keep_all = TRUE)

candidate_cols <- setdiff(colnames(clinical_patient), ignore_cols)

# Automatically classify clinical variables into Numeric, Binary, or Categorical
var_types <- map_dfr(candidate_cols, function(col_name) {
  raw_col <- clinical_patient[[col_name]]
  if (is.list(raw_col)) {
    raw_col <- map_chr(raw_col, ~ if (length(.x) == 0 || is.null(.x)) NA_character_ else as.character(.x[1]))
  }
  clean_data <- raw_col[!is.na(raw_col) & raw_col != "" & raw_col != "[Not Available]" & raw_col != "NA"]
  n_unique   <- length(unique(clean_data))
  
  if (n_unique < 2) return(tibble(Variable = col_name, Type = "Exclude"))
  is_pure_num <- is.numeric(raw_col) || (suppressWarnings(!any(is.na(as.numeric(clean_data)))) && n_unique > 15)
  
  type <- if (is_pure_num) "Numeric" else if (n_unique == 2) "Binary" else "Categorical"
  tibble(Variable = col_name, Type = type)
})

num_vars <- var_types %>% filter(Type == "Numeric") %>% pull(Variable)
bin_vars <- var_types %>% filter(Type == "Binary") %>% pull(Variable)
cat_vars <- var_types %>% filter(Type == "Categorical") %>% pull(Variable)

all_factors <- paste0("Factor", 1:n_factors)

all_factors_results <- map_dfr(all_factors, function(f_name) {
  curr_f <- factors_long %>%
    filter(Factor == f_name) %>%
    inner_join(clinical_patient, by = "patient_key")
  
  # A. Continuous covariates (Spearman correlation)
  res_num <- map_dfr(num_vars, function(var) {
    sub_df <- tibble(val = curr_f$factor_score, x = suppressWarnings(as.numeric(curr_f[[var]]))) %>%
      filter(!is.na(x), !is.na(val))
    if (nrow(sub_df) < 10) return(NULL)
    t <- cor.test(sub_df$val, sub_df$x, method = "spearman", exact = FALSE)
    tibble(Factor = f_name, Variable = var, Test = "Spearman", Stat_Type = "Rho", 
           Statistic = unname(t$estimate), p_value = t$p.value)
  })
  
  # B. Binary covariates (Wilcoxon rank-sum test)
  res_bin <- map_dfr(bin_vars, function(var) {
    sub_df <- tibble(val = curr_f$factor_score, grp = as.character(curr_f[[var]])) %>%
      filter(!is.na(grp), grp != "", grp != "[Not Available]", grp != "NA", !is.na(val)) %>%
      mutate(grp = droplevels(as.factor(grp)))
    if (nlevels(sub_df$grp) != 2 \vert{}\vert{} min(table(sub_df$grp)) < 2) return(NULL)
    t <- wilcox.test(val ~ grp, data = sub_df)
    tibble(Factor = f_name, Variable = var, Test = "Wilcoxon", Stat_Type = "W", 
           Statistic = unname(t$statistic), p_value = t$p.value)
  })
  
  # C. Multi-group categorical covariates (Kruskal-Wallis test)
  res_cat <- map_dfr(cat_vars, function(var) {
    sub_df <- tibble(val = curr_f$factor_score, grp = as.character(curr_f[[var]])) %>%
      filter(!is.na(grp), grp != "", grp != "[Not Available]", grp != "NA", !is.na(val))
    grp_counts   <- table(sub_df$grp)
    valid_groups <- names(grp_counts[grp_counts >= 2])
    sub_df <- sub_df %>% filter(grp %in% valid_groups) %>% mutate(grp = droplevels(as.factor(grp)))
    if (nlevels(sub_df$grp) < 2) return(NULL)
    t <- kruskal.test(val ~ grp, data = sub_df)
    tibble(Factor = f_name, Variable = var, Test = "Kruskal-Wallis", Stat_Type = "Chi-sq", 
           Statistic = unname(t$statistic), p_value = t$p.value)
  })
  
  bind_rows(res_num, res_bin, res_cat)
})

all_factors_tests <- all_factors_results %>%
  mutate(FDR_q = p.adjust(p_value, method = "BH"), Statistic = round(Statistic, 3)) %>%
  arrange(FDR_q, p_value)

write.csv(all_factors_tests, file.path(tab_dir, "mofa_clinical_associations_all.csv"), row.names = FALSE)
write.csv(filter(all_factors_tests, FDR_q < 0.05), file.path(tab_dir, "mofa_clinical_associations_sig.csv"), row.names = FALSE)
message("  -> Clinical association tables saved.")

message("[4/5] Characterizing Factor 1: Violin Plots & ROC Curves...")
# ------------------------------------------------------------------------------
f1_df <- data.frame(
  sample  = rownames(get_factors(model_final, factors = 1)[[1]]),
  Factor1 = get_factors(model_final, factors = 1)[[1]][, 1]
) %>%
  left_join(samples_metadata(model_final), by = "sample")

# A. Violin Plot by PAM50 Subtype
pam50_colors <- c("Basal" = "#E41A1C", "Her2" = "#FF7F00", "LumA" = "#33A02C", "LumB" = "#1F78B4", "Normal" = "pink")
f1_pam50_df  <- f1_df %>% filter(PAM50 %in% names(pam50_colors))

p_violin_pam50 <- ggplot(f1_pam50_df, aes(x = PAM50, y = Factor1, fill = PAM50)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = 21, outlier.size = 1.5) +
  scale_fill_manual(values = pam50_colors) +
  theme_classic(base_size = 12) +
  labs(x = "PAM50 Subtype", y = "Factor 1") +
  theme(legend.position = "none", axis.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", color = "black"))

ggsave(file.path(fig_dir, "Violin_PAM50.png"), plot = p_violin_pam50, width = 6, height = 4.5, dpi = 300)

# B. Violin Plot by TNBC Status
tnbc_colors <- c("Non-TNBC" = "#4575b4", "TNBC" = "darkred")
f1_tnbc_df  <- f1_df %>%
  filter(!is.na(tnbc)) %>%
  mutate(tnbc = factor(tnbc, levels = c("Non-TNBC", "TNBC")))

p_violin_tnbc <- ggplot(f1_tnbc_df, aes(x = tnbc, y = Factor1, fill = tnbc)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.18, fill = "white", color = "black", outlier.shape = 21, outlier.size = 1.5) +
  scale_fill_manual(values = tnbc_colors) +
  theme_classic(base_size = 12) +
  labs(x = "TNBC Status", y = "Factor 1") +
  theme(legend.position = "none", axis.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", color = "black"))

ggsave(file.path(fig_dir, "Violin_TNBC.png"), plot = p_violin_tnbc, width = 5.5, height = 4.5, dpi = 300)

# C. ROC Curves (Factor 1 inverted since Basal/TNBC lies on the negative pole)
roc_basal_df <- f1_pam50_df %>%
  mutate(D
