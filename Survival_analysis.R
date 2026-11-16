#!/usr/bin/env Rscript
# ==============================================================================
# Script: Survival_analysis.R
# Description: Performs Overall Survival (OS) analysis on TCGA-BRCA:
#                1) Multivariable Cox proportional hazards model on all MOFA+
#                   latent factors + Schoenfeld proportional hazards diagnostics
#                2) Kaplan-Meier survival curves (multivariable risk score and
#                   median splits for Factors 1, 4, and 13)
#                3) Univariate Cox models on top 100 continuous features for
#                   Factors 1, 4, and 13 (with BH FDR adjustment)
#                4) Multivariable Cox models on top 15 features per factor
#                   (with VIF multicollinearity and Schoenfeld checks)
#                5) Penalized Cox regression via Elastic Net (glmnet, alpha=0.5)
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(survival)
  library(survminer)
  library(broom)
  library(glmnet)
  library(ggplot2)
})

# Define directories
data_dir  <- file.path("data")
model_dir <- file.path("results", "models")
fig_dir   <- file.path("results", "figures")
tab_dir   <- file.path("results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/6] Loading Optimal MOFA+ Model and Preparing Survival Data...")
# ------------------------------------------------------------------------------
model_path <- file.path(model_dir, "best_mofa_model.hdf5")
clin_path  <- file.path(data_dir, "clinical_master.rds")

if (!file.exists(model_path)) stop("Missing 'best_mofa_model.hdf5'. Run 04_mofa_train.R first.")
if (!file.exists(clin_path))  stop("Missing 'clinical_master.rds'. Run 03_data_visualization.R first.")

model_final     <- load_model(model_path)
clinical_master <- readRDS(clin_path)

# Clean feature suffixes if present
features_names(model_final) <- lapply(features_names(model_final), function(feats) {
  sub("_(mRNA|CNV|Methylation|lncRNA|miRNA|Mutations)$", "", feats)
})

# Construct patient-level Overall Survival (OS) endpoint (12-char barcode)
survival_table <- clinical_master %>%
  mutate(barcode = substr(sample, 1, 12)) %>%
  distinct(barcode, .keep_all = TRUE) %>%
  mutate(
    status = case_when(
      vital_status == "Dead"  ~ 1,
      vital_status == "Alive" ~ 0,
      TRUE                    ~ NA_real_
    ),
    time = as.numeric(coalesce(days_to_death, days_to_last_follow_up))
  ) %>%
  filter(!is.na(status), !is.na(time), time > 0) %>%
  dplyr::select(barcode, time, status)

message(sprintf("  -> Patients with valid OS data: %d (Events/Deaths = %d, Censored = %d)",
                nrow(survival_table),
                sum(survival_table$status == 1),
                sum(survival_table$status == 0)))

# Extract patient-level latent factor scores (averaging duplicate aliquots if any)
n_factors   <- get_dimensions(model_final)$K
factor_cols <- paste0("Factor", 1:n_factors)

factors_patient <- as.data.frame(get_factors(model_final, factors = 1:n_factors)[[1]]) %>%
  rownames_to_column("sample") %>%
  mutate(barcode = substr(sample, 1, 12)) %>%
  group_by(barcode) %>%
  summarise(across(all_of(factor_cols), mean, na.rm = TRUE), .groups = "drop")

cox_data <- survival_table %>%
  inner_join(factors_patient, by = "barcode")

message(sprintf("  -> Matched survival + MOFA+ cohort: %d patients across %d factors.",
                nrow(cox_data), n_factors))

message("[2/6] Fitting Multivariable Cox Model Across All Latent Factors...")
# ------------------------------------------------------------------------------
formula_factors <- as.formula(paste("Surv(time, status) ~", paste(factor_cols, collapse = " + ")))
res_cox_factors <- coxph(formula_factors, data = cox_data)

tidy_cox_factors <- tidy(res_cox_factors, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::select(Factor = term, HR = estimate, Lower_95CI = conf.low,
                Upper_95CI = conf.high, z_stat = statistic, p_value = p.value)

write.csv(tidy_cox_factors, file.path(tab_dir, "cox_multivariable_factors.csv"), row.names = FALSE)
print(as.data.frame(tidy_cox_factors))

# Proportional Hazards Assumption Check (Scaled Schoenfeld Residuals)
ph_test <- cox.zph(res_cox_factors)
ph_df   <- as.data.frame(ph_test$table) %>% rownames_to_column("Term")
write.csv(ph_df, file.path(tab_dir, "cox_factors_schoenfeld_ph_test.csv"), row.names = FALSE)
message(sprintf("  -> Global PH Schoenfeld test p-value: %.4f", ph_test$table["GLOBAL", "p"]))

# Forest Plot of Factor Hazard Ratios
p_forest <- ggplot(tidy_cox_factors, aes(x = HR, y = factor(Factor, levels = rev(factor_cols)))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_errorbarh(aes(xmin = Lower_95CI, xmax = Upper_95CI), height = 0.25, color = "grey40") +
  geom_point(aes(color = p_value < 0.05), size = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#2c7fb8"),
                     labels = c("p >= 0.05", "p < 0.05"), name = "Significance") +
  theme_classic(base_size = 12) +
  labs(x = "Hazard Ratio (95% CI, Log Scale)", y = "") +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "Cox_factors_forest.png"), plot = p_forest, width = 6.5, height = 5, dpi = 300)

message("[3/6] Generating Kaplan-Meier Curves (Risk Score & Factors 1, 4, 13)...")
# ------------------------------------------------------------------------------
# A. Stratification by Multivariable Cox Linear Predictor (Risk Score)
cox_data$risk_score <- predict(res_cox_factors, type = "lp")
cox_data$risk_group <- factor(
  ifelse(cox_data$risk_score >= median(cox_data$risk_score, na.rm = TRUE), "High Risk", "Low Risk"),
  levels = c("Low Risk", "High Risk")
)

save_km_plot <- function(df, group_col, title_str, legend_labels, palette_cols, filename) {
  form <- as.formula(paste("Surv(time, status) ~", group_col))
  fit  <- surv_fit(form, data = df)
  
  p_km <- ggsurvplot(
    fit,
    data        = df,
    pval        = TRUE,
    conf.int    = TRUE,
    risk.table  = FALSE,
    legend.labs = legend_labels,
    palette     = palette_cols,
    title       = title_str,
    xlab        = "Time (Days)",
    ylab        = "Overall Survival Probability",
    ggtheme     = theme_classic(base_size = 12)
  )
  ggsave(file.path(fig_dir, filename), plot = p_km$plot, width = 6, height = 5, dpi = 300)
}

save_km_plot(
  df            = cox_data,
  group_col     = "risk_group",
  title_str     = "Multivariable Cox Risk Score",
  legend_labels = c("Low Risk", "High Risk"),
  palette_cols  = c("#2E9FDF", "#E64B35"),
  filename      = "KM_Cox_risk_score.png"
)

# B. Median Splits for Factor 1 (Subtype Axis) and Prognostic Factors 4 & 13
target_km_factors <- intersect(c("Factor1", "Factor4", "Factor13"), factor_cols)
for (f_col in target_km_factors) {
  grp_col <- paste0(f_col, "_group")
  cox_data[[grp_col]] <- factor(
    ifelse(cox_data[[f_col]] >= median(cox_data[[f_col]], na.rm = TRUE), "High", "Low"),
    levels = c("Low", "High")
  )
  save_km_plot(
    df            = cox_data,
    group_col     = grp_col,
    title_str     = gsub("Factor", "Factor ", f_col),
    legend_labels = c(paste("Low", gsub("Factor", "Factor ", f_col)),
                      paste("High", gsub("Factor", "Factor ", f_col))),
    palette_cols  = c("#E7B800", "#2E9FDF"),
    filename      = sprintf("KM_%s.png", tolower(f_col))
  )
}
message("  -> Kaplan-Meier plots saved to results/figures/.")

message("[4/6] Running Univariate Cox on Top 100 Continuous Features (Factors 1, 4, 13)...")
# ------------------------------------------------------------------------------
continuous_views <- c("mRNA", "lncRNA", "miRNA", "Methylation", "CNV")

# Extract patient-level continuous matrices
patient_matrices <- lapply(continuous_views, function(vw) {
  mat <- get_data(model_final, views = vw)[[1]][[1]]
  rownames(mat) <- sub(paste0("_", vw, "$"), "", rownames(mat))
  df <- as.data.frame(t(mat)) %>%
    rownames_to_column("sample") %>%
    mutate(barcode = substr(sample, 1, 12)) %>%
    group_by(barcode) %>%
    summarise(across(-sample, ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  return(df)
})
names(patient_matrices) <- continuous_views

# Helper to rank top features across all 5 continuous views by absolute factor weight
get_top_factor_features <- function(model, factor_idx, top_n = 100) {
  all_w <- map_dfr(continuous_views, function(vw) {
    w_vec <- get_weights(model, views = vw, factors = factor_idx)[[1]][, 1]
    tibble(
      View       = vw,
      Feature    = sub(paste0("_", vw, "$"), "", names(w_vec)),
      Weight     = as.numeric(w_vec),
      Abs_Weight = abs(as.numeric(w_vec))
    )
  })
  all_w %>% arrange(desc(Abs_Weight)) %>% slice_head(n = top_n)
}

run_univariate_top_features <- function(factor_idx, top_n = 100) {
  top_feats <- get_top_factor_features(model_final, factor_idx, top_n = top_n)
  
  res_uni <- pmap_dfr(top_feats, function(View, Feature, Weight, Abs_Weight) {
    view_df <- patient_matrices[[View]]
    if (!Feature %in% colnames(view_df)) return(NULL)
    
    sub_df <- survival_table %>%
      inner_join(view_df %>% dplyr::select(barcode, val = all_of(Feature)), by = "barcode") %>%
      filter(!is.na(val))
    
    if (nrow(sub_df) < 20 || var(sub_df$val) == 0) return(NULL)
    
    fit <- coxph(Surv(time, status) ~ val, data = sub_df)
    s   <- summary(fit)
    
    tibble(
      Factor     = paste0("Factor", factor_idx),
      View       = View,
      Feature    = Feature,
      Weight     = Weight,
      HR         = s$conf.int[1, "exp(coef)"],
      Lower_95CI = s$conf.int[1, "lower .95"],
      Upper_95CI = s$conf.int[1, "upper .95"],
      p_value    = s$coefficients[1, "Pr(>|z|)"]
    )
  })
  
  res_uni <- res_uni %>%
    mutate(FDR_q = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_q, p_value)
  
  write.csv(res_uni, file.path(tab_dir, sprintf("cox_univariate_Factor%d_top%d.csv", factor_idx, top_n)),
            row.names = FALSE)
  message(sprintf("  -> Factor %d Univariate Cox: %d/%d features significant at FDR < 0.05",
                  factor_idx, sum(res_uni$FDR_q < 0.05), nrow(res_uni)))
  return(res_uni)
}

target_factor_nums <- intersect(c(1, 4, 13), 1:n_factors)
uni_results_list   <- lapply(target_factor_nums, run_univariate_top_features, top_n = 100)

message("[5/6] Fitting Multivariable Cox Models on Top 15 Features (EPV >= 10)...")
# ------------------------------------------------------------------------------
# Helper to compute Variance Inflation Factors (VIF) for continuous predictors
compute_vif <- function(X_df) {
  cor_mat <- cor(X_df, use = "pairwise.complete.obs")
  inv_cor <- tryCatch(solve(cor_mat), error = function(e) NULL)
  if (is.null(inv_cor)) return(rep(NA_real_, ncol(X_df)))
  vif_vals <- diag(inv_cor)
  names(vif_vals) <- colnames(X_df)
  return(vif_vals)
}

run_multivariable_top15 <- function(factor_idx, top_n = 15) {
  top15 <- get_top_factor_features(model_final, factor_idx, top_n = top_n)
  
  # Assemble combined predictor matrix across views
  merged_df <- survival_table
  predictor_cols <- c()
  
  for (i in seq_len(nrow(top15))) {
    vw   <- top15$View[i]
    feat <- top15$Feature[i]
    col_id <- make.unique(c(predictor_cols, paste0(feat, "_", vw)))[length(predictor_cols) + 1]
    
    view_sub <- patient_matrices[[vw]] %>%
      dplyr::select(barcode, !!col_id := all_of(feat))
    
    merged_df <- inner_join(merged_df, view_sub, by = "barcode")
    predictor_cols <- c(predictor_cols, col_id)
  }
  
  merged_df <- merged_df %>% drop_na(all_of(predictor_cols))
  
  # Check Multicollinearity (VIF) and Proportional Hazards
  vif_vals <- compute_vif(merged_df[, predictor_cols, drop = FALSE])
  max_vif  <- max(vif_vals, na.rm = TRUE)
  
  form_multi <- as.formula(paste("Surv(time, status) ~",
                                 paste(sprintf("`%s`", predictor_cols), collapse = " + ")))
  fit_multi  <- coxph(form_multi, data = merged_df)
  ph_multi   <- cox.zph(fit_multi)
  
  tidy_multi <- tidy(fit_multi, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      Clean_Term = gsub("`", "", term),
      View       = top15$View,
      Feature    = top15$Feature,
      VIF        = vif_vals[Clean_Term]
    ) %>%
    dplyr::select(View, Feature, HR = estimate, Lower_95CI = conf.low,
                  Upper_95CI = conf.high, z_stat = statistic, p_value = p.value, VIF) %>%
    arrange(HR)
  
  write.csv(tidy_multi, file.path(tab_dir, sprintf("cox_multivariable_Factor%d_top%d.csv", factor_idx, top_n)),
            row.names = FALSE)
  message(sprintf("  -> Factor %d Top-15 Multivariable Cox: Max VIF = %.2f | Global PH p = %.3f | Sig features = %d",
                  factor_idx, max_vif, ph_multi$table["GLOBAL", "p"], sum(tidy_multi$p_value < 0.05)))
}

lapply(target_factor_nums, run_multivariable_top15, top_n = 15)

message("[6/6] Running Penalized Cox Regression via Elastic Net (glmnet, alpha = 0.5)...")
# ------------------------------------------------------------------------------
# Combine all continuous views (mRNA, lncRNA, miRNA, Methylation, CNV)
combined_blocks <- lapply(continuous_views, function(vw) {
  df <- patient_matrices[[vw]]
  # Prefix feature names with view to prevent cross-view symbol collisions
  colnames(df)[-1] <- paste0(colnames(df)[-1], "__", vw)
  return(df)
})

full_omics_df <- Reduce(function(x, y) full_join(x, y, by = "barcode"), combined_blocks)

enet_merged <- survival_table %>%
  inner_join(full_omics_df, by = "barcode")

X_raw <- as.matrix(enet_merged %>% dplyr::select(-barcode, -time, -status))
y_surv <- Surv(time = enet_merged$time, event = enet_merged$status)

# 1. Remove features with > 20% missing values across patients
na_prop <- colMeans(is.na(X_raw))
X_filt  <- X_raw[, na_prop <= 0.20, drop = FALSE]

# 2. Mean-impute remaining missing values
for (j in seq_len(ncol(X_filt))) {
  nas <- is.na(X_filt[, j])
  if (any(nas)) X_filt[nas, j] <- mean(X_filt[, j], na.rm = TRUE)
}

# 3. Remove zero-variance features and standardize to zero mean / unit variance
col_vars <- apply(X_filt, 2, var)
X_filt   <- X_filt[, col_vars > 0, drop = FALSE]
X_scaled <- scale(X_filt)

message(sprintf("  -> Elastic Net input matrix: %d patients x %d standardized continuous features.",
                nrow(X_scaled), ncol(X_scaled)))

# 4. 10-Fold Cross-Validation optimizing Harrell's C-index (alpha = 0.5)
set.seed(42)
cv_enet <- cv.glmnet(
  x            = X_scaled,
  y            = y_surv,
  family       = "cox",
  alpha        = 0.5,
  nfolds       = 10,
  type.measure = "C"
)

# Save CV curve plot
png(file.path(fig_dir, "glmnet_C_index.png"), width = 1800, height = 1500, res = 300)
plot(cv_enet)
dev.off()

# Extract non-zero features at lambda.1se
coef_1se <- coef(cv_enet, s = "lambda.1se")
nz_idx   <- which(as.numeric(coef_1se) != 0)

enet_features_df <- tibble(
  Raw_ID      = rownames(coef_1se)[nz_idx],
  Coefficient = as.numeric(coef_1se)[nz_idx],
  HR          = exp(as.numeric(coef_1se)[nz_idx])
) %>%
  separate(Raw_ID, into = c("Feature", "View"), sep = "__", extra = "merge") %>%
  mutate(
    Direction = ifelse(Coefficient > 0, "Increased Risk (HR > 1)", "Protective (HR < 1)")
  ) %>%
  arrange(View, HR)

write.csv(enet_features_df, file.path(tab_dir, "glmnet_cox_selected_features_1se.csv"), row.names = FALSE)

c_min <- cv_enet$cvm[cv_enet$index["min", 1]]
c_1se <- cv_enet$cvm[cv_enet$index["1se", 1]]
message(sprintf("  -> Elastic Net lambda.min = %.4f (C-index = %.3f, %d features)",
                cv_enet$lambda.min, c_min, sum(as.numeric(coef(cv_enet, s = "lambda.min")) != 0)))
message(sprintf("  -> Elastic Net lambda.1se = %.4f (C-index = %.3f, %d features)",
                cv_enet$lambda.1se, c_1se, nrow(enet_features_df)))

# Plot Diverging Bar Chart of Selected Features at lambda.1se
if (nrow(enet_features_df) > 0) {
  enet_plot_df <- enet_features_df %>%
    mutate(Label = paste0(Feature, " (", View, ")")) %>%
    arrange(View, Coefficient) %>%
    mutate(Label = factor(Label, levels = unique(Label)))
  
  p_enet <- ggplot(enet_plot_df, aes(x = Coefficient, y = Label, fill = Direction)) +
    geom_col(width = 0.75) +
    facet_grid(View ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("Increased Risk (HR > 1)" = "#E64B35",
                                 "Protective (HR < 1)"     = "#4DBBD5")) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    theme_classic(base_size = 10) +
    labs(
      title    = sprintf("Elastic Net Cox Signature (lambda.1se = %.4f)", cv_enet$lambda.1se),
      x        = expression("Regularized Coefficient ("*beta*")"),
      y        = "",
      fill     = "Prognostic Effect"
    ) +
    theme(strip.text.y = element_text(angle = 0, face = "bold"),
          legend.position = "top")
  
  ggsave(file.path(fig_dir, "glmnet_features.png"), plot = p_enet, width = 8, height = 11, dpi = 300)
}

message("\nDone! All survival analysis tables and figures saved to results/.")
