#!/usr/bin/env Rscript
# ==============================================================================
# Script: Mofa_train.R
# Description: Configures Python environment for mofapy2, trains 10 independent
#              MOFA+ models across distinct random seeds, evaluates convergence
#              via ELBO, and exports the optimal model for downstream analysis.
# Author: Miguel Camacho
# ==============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(MOFA2)
})

# Define directories
data_dir  <- file.path("data")
model_dir <- file.path("results", "models")
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/4] Configuring Python / mofapy2 Environment...")
# ------------------------------------------------------------------------------
# Note: MOFA2 relies on the Python backend 'mofapy2'. By default, MOFA2 can use
# basilisk (use_basilisk = TRUE). If using a custom virtualenv, it is created
# once if missing:
env_name <- "mofa_final_env"
use_custom_venv <- FALSE  # Set to TRUE if you prefer manual virtualenv over basilisk

if (use_custom_venv) {
  if (!reticulate::virtualenv_exists(env_name)) {
    message(sprintf("  -> Creating virtual environment '%s'...", env_name))
    reticulate::install_python(version = "3.10")
    reticulate::virtualenv_create(env_name, python = "3.10")
    reticulate::virtualenv_install(env_name, packages = c("mofapy2"))
  }
  reticulate::use_virtualenv(env_name, required = TRUE)
}

message("[2/4] Loading Aligned Multi-Omics Matrices...")
# ------------------------------------------------------------------------------
input_path <- file.path(data_dir, "preprocessed_mofa_inputs.rds")
if (!file.exists(input_path)) {
  stop("Missing 'preprocessed_mofa_inputs.rds' in data/. Run 01_preprocessing.R first.")
}

mofa_data <- readRDS(input_path)
mofa_data <- lapply(mofa_data, as.matrix)

# Verify view order and dimensions
message("  -> Loaded views: ", paste(names(mofa_data), collapse = ", "))
print(lapply(mofa_data, dim))

message("[3/4] Training 10 Independent MOFA+ Models Across Random Seeds...")
# ------------------------------------------------------------------------------
# Define 10 reproducible seeds (including seed 7)
seeds <- c(7, 14, 21, 42, 77, 100, 123, 256, 512, 1024)
trained_models <- list()
elbo_values    <- numeric(length(seeds))

for (i in seq_along(seeds)) {
  current_seed <- seeds[i]
  message(sprintf("\n--- Training Model %d/%d (Seed = %d) ---", i, length(seeds), current_seed))
  
  # Initialize fresh MOFA object
  MOFAobject <- create_mofa(mofa_data)
  
  # Data options: scale views to unit variance
  data_opts <- get_default_data_options(MOFAobject)
  data_opts$scale_views <- TRUE
  
  # Model options: K = 20 factors, two-level sparsity (ARD + Spike-and-Slab)
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors       <- 20
  model_opts$spikeslab_weights <- TRUE
  model_opts$ard_weights       <- TRUE
  model_opts$spikeslab_factors <- FALSE
  model_opts$ard_factors       <- FALSE
  model_opts$likelihoods       <- c("gaussian", "gaussian", "gaussian", 
                                    "gaussian", "gaussian", "bernoulli")
  
  # Training options: slow convergence (delta ELBO < 0.00005%), drop < 2% variance
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$drop_factor_threshold <- 0.02
  train_opts$convergence_mode      <- "slow"
  train_opts$maxiter               <- 2000
  train_opts$seed                  <- current_seed
  
  # Prepare MOFA object
  MOFAobject <- prepare_mofa(
    object           = MOFAobject,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )
  
  # Run MOFA+ and save seed-specific model
  outfile <- file.path(model_dir, sprintf("mofa_model_seed_%d.hdf5", current_seed))
  model_trained <- run_mofa(MOFAobject, outfile = outfile, use_basilisk = !use_custom_venv)
  
  trained_models[[i]] <- model_trained
  elbo_values[i]    <- get_elbo(model_trained)
  message(sprintf("  -> Seed %d converged with ELBO = %.2f", current_seed, elbo_values[i]))
}

message("\n[4/4] Selecting Optimal Model Based on Highest ELBO...")
# ------------------------------------------------------------------------------
summary_df <- data.frame(
  Run  = seq_along(seeds),
  Seed = seeds,
  ELBO = elbo_values
)
print(summary_df)

best_idx   <- which.max(elbo_values)
best_seed  <- seeds[best_idx]
best_file  <- file.path(model_dir, sprintf("mofa_model_seed_%d.hdf5", best_seed))
final_file <- file.path(model_dir, "best_mofa_model.hdf5")

file.copy(from = best_file, to = final_file, overwrite = TRUE)
write.csv(summary_df, file.path(model_dir, "mofa_seed_elbo_comparison.csv"), row.names = FALSE)

message(sprintf("Optimal model (Seed = %d, ELBO = %.2f) saved to: %s", 
                best_seed, elbo_values[best_idx], final_file))
message("Done!")
