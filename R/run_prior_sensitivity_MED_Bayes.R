# -----------------------------------------------------------------------------
# Targeted prior-sensitivity analysis for the conditional medoid Bayesian workflow
#
# Purpose:
#   Refit the decisive M2 vs M5 contrast on the already selected CART-BB medoid
#   completions under alternative regression-prior variances tau2_beta.
#
# Does NOT rerun the full master workflow and does NOT regenerate simulated data.
# -----------------------------------------------------------------------------

# Required objects before sourcing this file:
#   FULL_OUT           e.g. "outputs/full100_20260623_125642"
#   SOURCE_INDEX_FILE  e.g. "data/simulated_full100_20260623_125642/simulation_index.csv"

if (!exists("FULL_OUT")) {
  full_dirs <- list.dirs("outputs", recursive = FALSE, full.names = TRUE)
  FULL_OUT <- tail(sort(full_dirs[grepl("^full100_", basename(full_dirs))]), 1)
}

if (!exists("SOURCE_INDEX_FILE")) {
  idx_candidates <- list.files("data", pattern = "simulation_index\\.csv$", recursive = TRUE, full.names = TRUE)
  idx_candidates <- idx_candidates[grepl("simulated_full100", idx_candidates)]
  SOURCE_INDEX_FILE <- tail(sort(idx_candidates), 1)
}

message("FULL_OUT = ", FULL_OUT)
message("SOURCE_INDEX_FILE = ", SOURCE_INDEX_FILE)

stopifnot(dir.exists(FULL_OUT))
stopifnot(file.exists(SOURCE_INDEX_FILE))

source("R/07_julia_bridge.R")
source("R/17_sensitivity_summaries.R")

# User switches ----------------------------------------------------------------
if (!exists("RUN_SMOKE")) RUN_SMOKE <- FALSE
# Baseline tau2_beta = 1.0 is already in the main full100 results. Running it
# again is useful for a homogeneous sensitivity table but costs time.
if (!exists("RUN_BASELINE_TAU")) RUN_BASELINE_TAU <- FALSE

prior_grid <- if (RUN_BASELINE_TAU) c(0.25, 1.00, 4.00) else c(0.25, 4.00)
models_keep <- c("M2_true", "M5_full")
mechanisms_keep <- c("MCAR", "MAR", "MNAR")
reps_keep <- if (RUN_SMOKE) 1L else 1:100

# Use retained medoid completions as fixed-design input -------------------------
idx <- utils::read.csv(SOURCE_INDEX_FILE, stringsAsFactors = FALSE)
rep_col <- if ("rep_id" %in% names(idx)) "rep_id" else if ("rep" %in% names(idx)) "rep" else stop("No rep_id/rep column in SOURCE_INDEX_FILE.")
mech_col <- if ("mechanism" %in% names(idx)) "mechanism" else stop("No mechanism column in SOURCE_INDEX_FILE.")

idx[[mech_col]] <- toupper(as.character(idx[[mech_col]]))
idx <- idx[idx[[rep_col]] %in% reps_keep & idx[[mech_col]] %in% mechanisms_keep, , drop = FALSE]

medoid_dir <- file.path(FULL_OUT, "99_raw", "julia_cartbb", "medoid_imputations")
if (!dir.exists(medoid_dir)) stop("Medoid directory not found: ", medoid_dir)

medoid_path <- file.path(
  medoid_dir,
  sprintf("rep_%s_%s_medoid.csv", idx[[rep_col]], idx[[mech_col]])
)

missing_medoid <- medoid_path[!file.exists(medoid_path)]
if (length(missing_medoid) > 0L) {
  stop("Missing medoid files, first examples:\n", paste(head(missing_medoid, 10), collapse = "\n"))
}

idx$file_mis <- medoid_path
idx$file_full <- medoid_path
idx$file_cc <- medoid_path

sens_dir <- file.path(FULL_OUT, "05_sensitivity")
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)

medoid_index_file <- file.path(sens_dir, "medoid_index_for_prior_sensitivity.csv")
utils::write.csv(idx, medoid_index_file, row.names = FALSE)

model_file <- file.path(dirname(SOURCE_INDEX_FILE), "model_specs.csv")
if (!file.exists(model_file)) model_file <- "model_specs.csv"
stopifnot(file.exists(model_file))

# MCMC settings: keep this aligned with the main run if time permits.
# For a fast smoke test, these are overwritten below.
n_iter <- if (RUN_SMOKE) 1500L else 12000L
burnin <- if (RUN_SMOKE) 500L else 3000L
thin <- if (RUN_SMOKE) 10L else 6L

all_prior_results <- list()

for (tau2 in prior_grid) {
  tau_lab <- gsub("\\.", "p", format(tau2, trim = TRUE, scientific = FALSE))
  out_dir <- file.path(sens_dir, paste0("prior_tau2_", tau_lab))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("\nRunning prior sensitivity tau2_beta = ", tau2)

  run_julia_bayes(
    index_file = medoid_index_file,
    model_file = model_file,
    julia_project = "julia",
    julia_script = file.path("julia", "run_bayes.jl"),
    out_dir = out_dir,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    tau2_beta = tau2,
    c0 = 2.0,
    d0 = 1.0,
    impute_mode = "none",
    data_version = "mis",
    verbose_every = if (RUN_SMOKE) 500L else 2000L,
    models = models_keep,
    reps = reps_keep,
    mechanisms_filter = mechanisms_keep,
    expected_n_rep = if (RUN_SMOKE) 1L else 100L,
    expected_mechanisms = mechanisms_keep,
    export_draws = FALSE,
    diagnostics_dir = file.path(FULL_OUT, "06_mcmc_diagnostics")
  )

  f <- file.path(out_dir, "bayes_model_results.csv")
  if (!file.exists(f)) stop("Expected Julia result not found: ", f)

  res <- readr::read_csv(f, show_col_types = FALSE)
  if (!"workflow" %in% names(res)) res$workflow <- res$method
  res$workflow <- "MED-Bayes prior sensitivity"
  res$prior_tau2 <- tau2
  res$sensitivity_type <- "prior_tau2_beta"

  all_prior_results[[as.character(tau2)]] <- res
}

prior_evidence <- dplyr::bind_rows(all_prior_results)

prior_file <- file.path(sens_dir, "prior_evidence_by_replication.csv")
readr::write_csv(prior_evidence, prior_file)
message("Wrote: ", prior_file)

# Create the summary files expected by the paper-table scripts.
make_prior_sensitivity_tau2(FULL_OUT, prior_file = prior_file)

message("Prior-sensitivity run finished.")
