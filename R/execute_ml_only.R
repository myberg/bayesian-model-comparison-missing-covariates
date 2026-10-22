# -----------------------------------------------------------------------------
# execute_ml_only.R
# Reruns only the ML benchmarks on the exact CSV inputs exported for a
# completed Bayesian run. The Julia/Bayesian workflow is never called.
# -----------------------------------------------------------------------------

# Preserve user supplied settings before clearing the workspace.
# You can set them either by editing the two objects below or by calling
# Sys.setenv(RUN_MODE = "full100", SOURCE_INDEX_FILE = "...") before sourcing.
.user_run_mode <- Sys.getenv("RUN_MODE", unset = "")
.user_source_index_file <- Sys.getenv("SOURCE_INDEX_FILE", unset = "")
rm(list = setdiff(ls(), c(".user_run_mode", ".user_source_index_file")))

source("R/00_config.R")
load_required_packages(c(required_packages, "missForest"))
source("R/01_simulation_dgp_missingness.R")
source("R/04_ml_benchmarks.R")
source("R/05_run_simulation.R")
source("R/08_ml_before_deletion.R")
source("R/09_ml_imputation_estimates.R")
source("R/10_mice_d2_lrt.R")

RUN_MODE <- if (nzchar(.user_run_mode)) .user_run_mode else "full100"
stopifnot(RUN_MODE %in% c("pilot10", "full50a", "full50b", "full100"))

# Set this to the simulation_index.csv belonging to the completed Gibbs run.
# You may either edit this fallback value here or set it before sourcing via:
# Sys.setenv(SOURCE_INDEX_FILE = "data/simulated_full100_.../simulation_index.csv")
SOURCE_INDEX_FILE <- if (nzchar(.user_source_index_file)) {
  .user_source_index_file
} else {
  "data/simulated_full100_20260623_125642/simulation_index.csv"
}

if (!nzchar(SOURCE_INDEX_FILE) || !file.exists(SOURCE_INDEX_FILE)) {
  candidates <- list.files(
    "data", pattern = "^simulation_index[.]csv$",
    recursive = TRUE, full.names = TRUE
  )
  stop(
    "Set SOURCE_INDEX_FILE to the simulation_index.csv of the completed Gibbs run.\n",
    "Available index files:\n", paste(candidates, collapse = "\n"),
    call. = FALSE
  )
}

rep_ids <- switch(
  RUN_MODE,
  pilot10 = seq_len(10L),
  full50a = seq_len(50L),
  full50b = 51:100,
  full100 = seq_len(100L)
)

mechanisms <- c("MCAR", "MAR", "MNAR")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_root <- file.path("outputs", paste0("ml_only_", RUN_MODE, "_", timestamp))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

idx <- utils::read.csv(SOURCE_INDEX_FILE, stringsAsFactors = FALSE)
required_index_columns <- c("rep_id", "mechanism", "file_full", "file_mis")
missing_index_columns <- setdiff(required_index_columns, names(idx))
if (length(missing_index_columns) > 0L) {
  stop("The simulation index is missing columns: ",
       paste(missing_index_columns, collapse = ", "), call. = FALSE)
}
idx <- idx[idx$rep_id %in% rep_ids & idx$mechanism %in% mechanisms, , drop = FALSE]
expected_rows <- length(rep_ids) * length(mechanisms)
if (nrow(idx) != expected_rows) {
  stop("Expected ", expected_rows, " index rows but found ", nrow(idx),
       ". Check RUN_MODE and SOURCE_INDEX_FILE.", call. = FALSE)
}

ml_rows <- vector("list", nrow(idx))
estimate_rows <- vector("list", nrow(idx))
d2_rows <- vector("list", nrow(idx))
runtime_rows <- vector("list", nrow(idx))

for (i in seq_len(nrow(idx))) {
  mis_file <- .resolve_index_path(idx$file_mis[i], SOURCE_INDEX_FILE)
  if (!file.exists(mis_file)) stop("Missing input file: ", mis_file, call. = FALSE)
  data_mis <- utils::read.csv(mis_file, stringsAsFactors = FALSE)
  if ("id" %in% names(data_mis)) data_mis$id <- factor(data_mis$id)

  one <- run_ml_benchmarks_one(
    sim_item = list(
      rep_id = as.integer(idx$rep_id[i]),
      mechanism = idx$mechanism[i],
      data_mis = data_mis
    ),
    model_specs = model_specs,
    run_missforest = TRUE,
    mice_m = 20L,
    mice_maxit = 10L,
    seed = 202605L
  )
  ml_rows[[i]] <- one$ml
  estimate_rows[[i]] <- one$ml_estimates
  d2_rows[[i]] <- one$mice_d2_tests
  runtime_rows[[i]] <- one$ml_runtime_components
}

res_ml <- list(
  ml = .bind_rows_base(ml_rows),
  ml_estimates = .bind_rows_base(estimate_rows),
  mice_d2_tests = .bind_rows_base(d2_rows),
  ml_runtime_components = .bind_rows_base(runtime_rows)
)

before <- run_ml_before_deletion_from_index(
  index_file = SOURCE_INDEX_FILE,
  model_specs_arg = model_specs,
  out_dir = file.path(out_root, "before_deletion"),
  models = names(model_specs),
  reps = rep_ids,
  mechanisms_filter = mechanisms
)
res_ml <- append_ml_before_deletion(res_ml, before)

utils::write.csv(res_ml$ml,
                 file.path(out_root, "ml_model_results.csv"), row.names = FALSE)
utils::write.csv(res_ml$ml_estimates,
                 file.path(out_root, "ml_estimates.csv"), row.names = FALSE)
utils::write.csv(res_ml$mice_d2_tests,
                 file.path(out_root, "mice_d2_nested_tests.csv"), row.names = FALSE)
utils::write.csv(res_ml$ml_runtime_components,
                 file.path(out_root, "ml_runtime_components.csv"), row.names = FALSE)

message("ML-only workflow finished: ", out_root)
