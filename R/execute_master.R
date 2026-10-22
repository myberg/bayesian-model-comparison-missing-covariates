# ---
# execute_master_v14.R
# Master workflow for the SMAP simulation runs.
#
# Purpose:
#   1. Run a smoke-test with one MCAR replication first.
#   2. Keep the Bayesian workflows in separate output directories so that
#      bayes_model_results.csv and related files cannot overwrite each other.
#   3. Export harmonised raw files into outputs/.
#   4. Prepare the central CSV structure used in Section 4 of the manuscript.
#
# IMPORTANT:
#   - Start with RUN_MODE <- "smoke1".
#   - Continue with RUN_MODE <- "pilot10" before the full C = 100 run.
#   - Use RUN_MODE <- "full50a" for replications 1--50.
#   - Use RUN_MODE <- "full50b" later for replications 51--100.
#   - RUN_MODE <- "full100" remains available for one uninterrupted run.
#   - If MNAR is not yet implemented in R/01_simulation_dgp_missingness.R,
#     temporarily set INCLUDE_MNAR <- FALSE.
# ---

# clear everything, start new
rm(list = ls())

# 0. Load project code---------------------------------------------------------

source("R/00_config.R")
load_required_packages()

source("R/01_simulation_dgp_missingness.R")
source("R/02_cart_bb_update.R")
source("R/03_bayes_sampler.R")
source("R/04_ml_benchmarks.R")
source("R/05_run_simulation.R")
source("R/06_collect_tables.R")
source("R/07_julia_bridge.R")
source("R/08_ml_before_deletion.R")
source("R/09_ml_imputation_estimates.R")
source("R/10_mice_d2_lrt.R")
source("R/11_export_bayes_complete_case.R")




# 1. Central setup-------------------------------------------------------------

RUN_MODE <- "full100"
stopifnot(RUN_MODE %in% c("smoke1", "pilot10", "full50a", "full50b", "full100"))

RUN_REP_IDS <- switch(
  RUN_MODE,
  smoke1 = 1L,
  pilot10 = seq_len(10L),
  full50a = seq_len(50L),
  full50b = 51:100,
  full100 = seq_len(100L)
)

SMOKE_TEST <- identical(RUN_MODE, "smoke1")
PILOT_TEST <- identical(RUN_MODE, "pilot10")
FULL_RUN <- RUN_MODE %in% c("full50a", "full50b", "full100")
INCLUDE_MNAR <- TRUE

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0(RUN_MODE, "_", timestamp)
run_tag <- run_id


rev_out <- file.path("outputs", run_id)
rev_root <- rev_out
raw_root <- file.path(rev_out, "99_raw")
data_root <- file.path("data", paste0("simulated_", run_id))

diagnostics_dir <- file.path(rev_out, "06_mcmc_diagnostics")

# setup diagnostics

export_mcmc_draws <- TRUE
diagnostic_reps <- min(RUN_REP_IDS)
diagnostic_models <- c("M2_true")
diagnostic_mechanisms <- if (SMOKE_TEST) {
  "MCAR"
} else if (INCLUDE_MNAR) {
  c("MCAR", "MAR", "MNAR")
} else {
  c("MCAR", "MAR")
}


dir.create(rev_out, recursive = TRUE, showWarnings = FALSE)


setup <- list(
  n_rep = length(RUN_REP_IDS),
  rep_ids = RUN_REP_IDS,
  mechanisms = if (SMOKE_TEST) {
    "MCAR"
  } else if (INCLUDE_MNAR) {
    c("MCAR", "MAR", "MNAR")
  } else {
    c("MCAR", "MAR")
  },

  run_ml = TRUE,
  run_missforest = TRUE,
  run_ml_before_deletion = TRUE,
  export_data_for_julia = TRUE,
  force_reexport = FALSE,

  run_julia_before_deletion = TRUE,
  run_julia_complete_case = TRUE,
  run_julia_cartbb = TRUE,

  mice_m = 20L,
  mice_maxit = 10L,

  julia_models = c(
    "M1_intercept", "M2_true", "M3_reduced",
    "M4_noise", "M5_full", "M6_alt"
  ),
  
  # julia_models = c("M2_true"),

  julia_n_iter = if (FULL_RUN) 20000L else switch(RUN_MODE, smoke1 = 1000L, pilot10 = 4000L),
  julia_burnin = if (FULL_RUN) 4000L else switch(RUN_MODE, smoke1 = 250L, pilot10 = 1000L),
  julia_thin = if (FULL_RUN) 2L else 5L,

  # julia_cart_every = 1L,
  # julia_cart_max_depth = 7L,
  # julia_cart_min_samples_leaf = 5L,
  # julia_cart_max_donors = 100L,
  
  julia_cart_every = 10L,
  julia_cart_max_depth = 5L,
  julia_cart_min_samples_leaf = 10L,
  julia_cart_max_donors = 50L,
  julia_medoid_reference_model = "M2_true",
  julia_medoid_thin = if (FULL_RUN) 20L else switch(RUN_MODE, smoke1 = 5L, pilot10 = 10L),
  julia_chib_reduced_iter = if (FULL_RUN) 4000L else switch(RUN_MODE, smoke1 = 500L, pilot10 = 2000L),
  julia_chib_reduced_burnin = if (FULL_RUN) 1000L else switch(RUN_MODE, smoke1 = 125L, pilot10 = 500L),
  julia_chib_reduced_thin = 2L,
  julia_verbose_every = if (FULL_RUN) 1000L else switch(RUN_MODE, smoke1 = 250L, pilot10 = 500L)
  
)

if (isTRUE(setup$run_missforest) && !requireNamespace("missForest", quietly = TRUE)) {
  stop("Install the optional package 'missForest' before running this workflow.",
       call. = FALSE)
}


# 2. Helper functions----------------------------------------------------------

make_dirs <- function(root) {
  dirs <- file.path(root, c(
    "00_metadata",
    "01_parameter_recovery",
    "02_model_selection",
    "03_evidence",
    "04_predictive_criteria",
    "05_sensitivity",
    "06_mcmc_diagnostics",
    "07_runtime",
    "08_figures_data",
    "99_raw/r_ml",
    "99_raw/julia_before_deletion",
    "99_raw/julia_complete_case",
    "99_raw/julia_cartbb"
  ))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(x) && is.data.frame(x) && nrow(x) > 0L) {
    utils::write.csv(x, file, row.names = FALSE)
  }
}

add_or_replace <- function(x, name, value) {
  if (!is.null(x) && is.data.frame(x) && nrow(x) > 0L) {
    x[[name]] <- value
  }
  x
}

tag_bayes_results <- function(res,
                              workflow,
                              data_version,
                              impute_mode,
                              evidence_type) {
  if (is.null(res)) return(res)

  for (nm in names(res)) {
    if (is.data.frame(res[[nm]]) && nrow(res[[nm]]) > 0L) {
      res[[nm]]$method <- workflow
      res[[nm]]$workflow <- workflow
      res[[nm]]$data_version <- data_version
      res[[nm]]$impute_mode <- impute_mode
      if (nm %in% c("bayes", "bayes_marginal_likelihood")) {
        res[[nm]]$evidence_type <- evidence_type
      }
    }
  }
  res
}

bind_existing <- function(x) {
  x <- Filter(function(z) !is.null(z) && is.data.frame(z) && nrow(z) > 0L, x)
  if (length(x) == 0L) return(NULL)
  do.call(rbind, x)
}

run_one_julia_workflow <- function(workflow,
                                   data_version,
                                   impute_mode,
                                   conditional_medoid = FALSE,
                                   result_data_version = if (conditional_medoid) "medoid" else data_version,
                                   result_impute_mode = if (conditional_medoid) "none" else impute_mode,
                                   out_dir,
                                   export_draws = export_draws,
                                   diagnostics_dir = diagnostics_dir,
                                   draw_reps = draw_reps,
                                   draw_models = draw_models,
                                   draw_mechanisms = draw_mechanisms,
                                   evidence_type,
                                   n_iter = setup$julia_n_iter,
                                   burnin = setup$julia_burnin,
                                   thin = setup$julia_thin) {
  message("Running Julia workflow: ", workflow)
  unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  run_julia_bayes(
    index_file = file.path(data_root, "simulation_index.csv"),
    model_file = file.path(data_root, "model_specs.csv"),
    julia_project = "julia",
    julia_script = file.path("julia", "run_bayes.jl"),
    out_dir = out_dir,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    impute_mode = impute_mode,
    data_version = data_version,
    cart_every = setup$julia_cart_every,
    cart_max_depth = setup$julia_cart_max_depth,
    cart_min_samples_leaf = setup$julia_cart_min_samples_leaf,
    cart_max_donors = setup$julia_cart_max_donors,
    conditional_medoid = conditional_medoid,
    medoid_reference_model = setup$julia_medoid_reference_model,
    medoid_n_iter = n_iter,
    medoid_burnin = burnin,
    medoid_thin = setup$julia_medoid_thin,
    chib_reduced_iter = setup$julia_chib_reduced_iter,
    chib_reduced_burnin = setup$julia_chib_reduced_burnin,
    chib_reduced_thin = setup$julia_chib_reduced_thin,
    verbose_every = setup$julia_verbose_every,
    models = setup$julia_models,
    reps = setup$rep_ids,
    expected_n_rep = setup$n_rep,
    expected_rep_ids = setup$rep_ids,
    expected_mechanisms = setup$mechanisms,
    
    export_draws = export_draws,
    diagnostics_dir = diagnostics_dir,
    draw_reps = draw_reps,
    draw_models = draw_models,
    draw_mechanisms = draw_mechanisms
  )

  res <- read_julia_bayes_results(out_dir)
  tag_bayes_results(
    res = res,
    workflow = workflow,
    data_version = result_data_version,
    impute_mode = result_impute_mode,
    evidence_type = evidence_type
  )
}

normalise_runtime <- function(x, source_label) {
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0L) return(NULL)
  if (!"runtime_seconds" %in% names(x) && "elapsed_seconds" %in% names(x)) {
    x$runtime_seconds <- x$elapsed_seconds
  }
  if (!"data_version" %in% names(x)) x$data_version <- NA_character_
  if (!"impute_mode" %in% names(x)) x$impute_mode <- NA_character_
  if (!"workflow" %in% names(x) && "method" %in% names(x)) x$workflow <- x$method
  x$source <- source_label
  keep <- intersect(
    c("source", "method", "workflow", "model", "rep", "rep_id", "mechanism",
      "data_version", "impute_mode", "runtime_seconds"),
    names(x)
  )
  x[, keep, drop = FALSE]
}

select_models_from_criterion <- function(df, criterion, direction = c("min", "max")) {
  direction <- match.arg(direction)
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L || !criterion %in% names(df)) {
    return(NULL)
  }
  rep_col <- if ("rep" %in% names(df)) "rep" else if ("rep_id" %in% names(df)) "rep_id" else NA_character_
  if (is.na(rep_col)) return(NULL)
  if (!"workflow" %in% names(df) && "method" %in% names(df)) df$workflow <- df$method
  if (!"mechanism" %in% names(df) || !"model" %in% names(df)) return(NULL)

  split_key <- interaction(df[[rep_col]], df$mechanism, df$workflow, drop = TRUE)
  pieces <- split(df, split_key)
  out <- lapply(pieces, function(z) {
    z <- z[!is.na(z[[criterion]]), , drop = FALSE]
    if (nrow(z) == 0L) return(NULL)
    ii <- if (direction == "min") which.min(z[[criterion]]) else which.max(z[[criterion]])
    data.frame(
      rep = z[[rep_col]][ii],
      mechanism = z$mechanism[ii],
      workflow = z$workflow[ii],
      criterion = criterion,
      selected_model = z$model[ii],
      stringsAsFactors = FALSE
    )
  })
  bind_existing(out)
}

selection_rates <- function(selected) {
  if (is.null(selected) || nrow(selected) == 0L) return(NULL)
  tab <- as.data.frame(table(
    selected$mechanism,
    selected$workflow,
    selected$criterion,
    selected$selected_model,
    useNA = "no"
  ), stringsAsFactors = FALSE)
  names(tab) <- c("mechanism", "workflow", "criterion", "selected_model", "n")
  tab <- tab[tab$n > 0L, , drop = FALSE]
  totals <- aggregate(n ~ mechanism + workflow + criterion, data = tab, sum)
  names(totals)[names(totals) == "n"] <- "n_total"
  tab <- merge(tab, totals, by = c("mechanism", "workflow", "criterion"), all.x = TRUE)
  tab$selection_rate <- tab$n / tab$n_total
  tab$mcse_selection_rate <- sqrt(
    tab$selection_rate * (1 - tab$selection_rate) / tab$n_total
  )
  tab
}

summarise_runtime <- function(runtime_long) {
  if (is.null(runtime_long) || nrow(runtime_long) == 0L) return(NULL)
  aggregate(
    runtime_seconds ~ source + workflow + mechanism + data_version + impute_mode,
    data = runtime_long,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q25 = unname(stats::quantile(x, 0.25, na.rm = TRUE)),
      q75 = unname(stats::quantile(x, 0.75, na.rm = TRUE)),
      sd = stats::sd(x, na.rm = TRUE)
    )
  ) |>
    do.call(what = data.frame)
}

# Metadata summaries. These are intentionally simple and robust.
write_metadata <- function() {
  simulation_settings <- data.frame(
    setting = c(
      "run_tag", "run_mode", "smoke_test", "n_rep", "rep_first", "rep_last", "mechanisms", "n_iter",
      "burnin", "thin", "mice_m", "mice_maxit", "cart_every", "cart_max_depth",
      "cart_min_samples_leaf", "cart_max_donors"
    ),
    value = c(
      run_tag, RUN_MODE, SMOKE_TEST, setup$n_rep, min(setup$rep_ids), max(setup$rep_ids),
      paste(setup$mechanisms, collapse = ";"),
      setup$julia_n_iter, setup$julia_burnin, setup$julia_thin,
      setup$mice_m, setup$mice_maxit,
      setup$julia_cart_every, setup$julia_cart_max_depth,
      setup$julia_cart_min_samples_leaf, setup$julia_cart_max_donors
    )
  )
  safe_write_csv(simulation_settings, file.path(rev_root, "00_metadata", "simulation_settings.csv"))

  workflow_dictionary <- data.frame(
    workflow = c(
      "Before deletion + ML",
      "Before deletion + Bayesian MCMC",
      "Complete case + ML",
      "Complete case + Bayesian MCMC",
      "MICE-CART + ML",
      "missForest + ML",
      "CART-BB medoid + conditional Bayesian MCMC"
    ),
    data_treatment = c(
      "Complete baseline data",
      "Complete baseline data",
      "Listwise deletion",
      "Listwise deletion",
      "Multiple imputation with MICE-CART",
      "Single imputation with missForest",
      "Common CART-BB medoid completion"
    ),
    estimator = c("ML", "Bayesian MCMC", "ML", "Bayesian MCMC", "ML", "ML", "Bayesian MCMC")
  )
  safe_write_csv(workflow_dictionary, file.path(rev_root, "00_metadata", "workflow_dictionary.csv"))

  if (exists("model_specs")) {
    candidate_models <- data.frame(
      model = names(model_specs),
      specification = vapply(model_specs, function(z) paste(deparse(z), collapse = " "), character(1)),
      stringsAsFactors = FALSE
    )
    safe_write_csv(candidate_models, file.path(rev_root, "00_metadata", "candidate_models.csv"))
  }

  software_versions <- data.frame(
    component = c("R", "platform"),
    version = c(as.character(getRversion()), R.version$platform),
    stringsAsFactors = FALSE
  )
  safe_write_csv(software_versions, file.path(rev_root, "00_metadata", "software_versions.csv"))
}



# 3. Create output structure and export simulation data------------------------

make_dirs(rev_root)
write_metadata()

if (isTRUE(setup$export_data_for_julia)) {
  if (isTRUE(setup$force_reexport)) {
    unlink(data_root, recursive = TRUE, force = TRUE)
  }

  sim_export <- export_simulation_data(
    n_rep = setup$n_rep,
    rep_ids = setup$rep_ids,
    mechanisms = setup$mechanisms,
    out_dir = data_root
  )

  inspect_simulation_index(
    sim_export$index_file,
    expected_n_rep = setup$n_rep,
    expected_rep_ids = setup$rep_ids,
    expected_mechanisms = setup$mechanisms,
    stop_on_mismatch = TRUE
  )

  file.copy(
    sim_export$index_file,
    file.path(rev_root, "00_metadata", "simulation_index.csv"),
    overwrite = TRUE
  )
}

idx <- utils::read.csv(file.path(data_root, "simulation_index.csv"))
print(table(idx$rep_id, idx$mechanism))






# 4. R-side ML, MICE-CART, missForest, and before-deletion ML------------------

if (isTRUE(setup$run_ml)) {
  res_ml <- run_simulation(
    n_rep = setup$n_rep,
    rep_ids = setup$rep_ids,
    mechanisms = setup$mechanisms,
    run_ml = TRUE,
    run_bayes = FALSE,
    run_missforest = setup$run_missforest,
    mice_m = setup$mice_m,
    mice_maxit = setup$mice_maxit
  )

  if (isTRUE(setup$run_ml_before_deletion)) {
    ml_before_deletion <- run_ml_before_deletion_from_index(
      index_file = file.path(data_root, "simulation_index.csv"),
      model_specs_arg = model_specs,
      out_dir = file.path(raw_root, "r_ml"),
      models = names(model_specs),
      reps = setup$rep_ids,
      mechanisms_filter = setup$mechanisms
    )
    res_ml <- append_ml_before_deletion(res_ml, ml_before_deletion)
  }

  safe_write_csv(res_ml$ml, file.path(raw_root, "r_ml", "ml_model_results.csv"))
  safe_write_csv(res_ml$ml_estimates, file.path(raw_root, "r_ml", "ml_estimates.csv"))
  safe_write_csv(res_ml$mice_d2_tests, file.path(raw_root, "r_ml", "mice_d2_nested_tests.csv"))
  safe_write_csv(res_ml$ml_runtime_components, file.path(raw_root, "r_ml", "ml_runtime_components.csv"))
} else {
  res_ml <- list(
    ml = NULL, ml_estimates = NULL,
    mice_d2_tests = NULL, ml_runtime_components = NULL
  )
}



# 5. Julia-side Bayesian workflows, each in its own output directory-----------

julia_before <- NULL
julia_cc <- NULL
julia_cartbb <- NULL


idx <- read.csv(file.path(data_root, "simulation_index.csv"))

actual_n_rep <- length(unique(idx$rep_id))
actual_mechanisms <- sort(unique(idx$mechanism))

message("Index contains ", actual_n_rep, " replications and mechanisms: ",
        paste(actual_mechanisms, collapse = ", "))

setup$n_rep <- actual_n_rep
setup$mechanisms <- actual_mechanisms


# before deletion
if (isTRUE(setup$run_julia_before_deletion)) {
  julia_before <- run_one_julia_workflow(
    workflow = "Before deletion + Bayesian MCMC",
    data_version = "full",
    impute_mode = "none",
    out_dir = file.path(raw_root, "julia_before_deletion"),
    evidence_type = "Fixed-design marginal likelihood via reduced-Gibbs Chib",
    export_draws = TRUE,
    diagnostics_dir = file.path(rev_out, "06_mcmc_diagnostics"),
    draw_reps = diagnostic_reps,
    draw_models = c("M2_true"),
    draw_mechanisms = diagnostic_mechanisms
  )
}

# complete cases
if (isTRUE(setup$run_julia_complete_case)) {
  export_common_complete_case_for_julia(
    index_file = file.path(data_root, "simulation_index.csv"),
    model_specs_arg = model_specs,
    id_var = "id",
    y_var = "y"
  )

  julia_cc <- run_one_julia_workflow(
    workflow = "Complete case + Bayesian MCMC",
    data_version = "cc",
    impute_mode = "none",
    out_dir = file.path(raw_root, "julia_complete_case"),
    evidence_type = "Fixed-design marginal likelihood via reduced-Gibbs Chib",
    export_draws = TRUE,
    diagnostics_dir = file.path(rev_out, "06_mcmc_diagnostics"),
    draw_reps = diagnostic_reps,
    draw_models = c("M2_true"),
    draw_mechanisms = diagnostic_mechanisms
  )
}


# CART-BB reference imputation, common medoid, then conditional evidence
if (isTRUE(setup$run_julia_cartbb)) {
  julia_cartbb <- run_one_julia_workflow(
    workflow = "CART-BB medoid + conditional Bayesian MCMC",
    data_version = "mis",
    impute_mode = "cartbb_every",
    conditional_medoid = TRUE,
    out_dir = file.path(raw_root, "julia_cartbb"),
    evidence_type = "Conditional marginal likelihood p(y|X_medoid,M) via reduced-Gibbs Chib",
    export_draws = TRUE,
    diagnostics_dir = file.path(rev_out, "06_mcmc_diagnostics"),
    draw_reps = diagnostic_reps,
    draw_models = c("M2_true"),
    draw_mechanisms = diagnostic_mechanisms
  )
}



# 6. Combine Bayesian result objects and write harmonized raw files------------
bayes_model_results <- bind_existing(list(
  if (!is.null(julia_before)) julia_before$bayes else NULL,
  if (!is.null(julia_cc)) julia_cc$bayes else NULL,
  if (!is.null(julia_cartbb)) julia_cartbb$bayes else NULL
))

bayes_estimates <- bind_existing(list(
  if (!is.null(julia_before)) julia_before$bayes_estimates else NULL,
  if (!is.null(julia_cc)) julia_cc$bayes_estimates else NULL,
  if (!is.null(julia_cartbb)) julia_cartbb$bayes_estimates else NULL
))

bayes_evidence <- bind_existing(list(
  if (!is.null(julia_before)) julia_before$bayes_marginal_likelihood else NULL,
  if (!is.null(julia_cc)) julia_cc$bayes_marginal_likelihood else NULL,
  if (!is.null(julia_cartbb)) julia_cartbb$bayes_marginal_likelihood else NULL
))

bayes_pred <- bind_existing(list(
  if (!is.null(julia_before)) julia_before$bayes_predictive_criteria else NULL,
  if (!is.null(julia_cc)) julia_cc$bayes_predictive_criteria else NULL,
  if (!is.null(julia_cartbb)) julia_cartbb$bayes_predictive_criteria else NULL
))

medoid_metadata <- if (!is.null(julia_cartbb)) {
  julia_cartbb$medoid_imputation_metadata
} else {
  NULL
}

safe_write_csv(bayes_model_results, file.path(rev_root, "02_model_selection", "bayes_model_results_all_workflows.csv"))
safe_write_csv(bayes_estimates, file.path(rev_root, "01_parameter_recovery", "bayes_estimates_all_workflows.csv"))
safe_write_csv(bayes_evidence, file.path(rev_root, "03_evidence", "evidence_by_replication.csv"))
safe_write_csv(bayes_pred, file.path(rev_root, "04_predictive_criteria", "predictive_criteria_by_replication.csv"))
safe_write_csv(medoid_metadata, file.path(rev_root, "03_evidence", "medoid_imputation_metadata.csv"))


# Also keep the ML raw files in their target folders.
safe_write_csv(res_ml$ml, file.path(rev_root, "02_model_selection", "ml_model_results_all_workflows.csv"))
safe_write_csv(res_ml$ml_estimates, file.path(rev_root, "01_parameter_recovery", "ml_estimates_all_workflows.csv"))
safe_write_csv(res_ml$mice_d2_tests, file.path(rev_root, "02_model_selection", "mice_d2_nested_tests.csv"))
safe_write_csv(res_ml$ml_runtime_components, file.path(rev_root, "07_runtime", "ml_runtime_components.csv"))



# 7. Model-selection and runtime summaries-------------------------------------

selected_models <- bind_existing(list(
  select_models_from_criterion(res_ml$ml, "BIC", direction = "min"),
  select_models_from_criterion(res_ml$ml, "AIC", direction = "min"),
  select_models_from_criterion(bayes_model_results, "log_marginal_likelihood", direction = "max"),
  select_models_from_criterion(bayes_model_results, "WAIC", direction = "min"),
  select_models_from_criterion(bayes_model_results, "LOO", direction = "min")
))


safe_write_csv(selected_models, file.path(rev_root, "02_model_selection", "selected_models_by_replication.csv"))
safe_write_csv(selection_rates(selected_models), file.path(rev_root, "02_model_selection", "model_selection_rates_all_mechanisms.csv"))

runtime_long <- bind_existing(list(
  normalise_runtime(res_ml$ml, "ML"),
  normalise_runtime(bayes_model_results, "Bayes")
))

safe_write_csv(runtime_long, file.path(rev_root, "07_runtime", "runtime_by_model.csv"))
safe_write_csv(summarise_runtime(runtime_long), file.path(rev_root, "07_runtime", "runtime_by_workflow_mechanism.csv"))

runtime_context <- data.frame(
  language = c("R", "Julia"),
  implementation = c("ML/MICE/missForest", "Bayesian MCMC"),
  n_iter = c(NA, setup$julia_n_iter),
  burnin = c(NA, setup$julia_burnin),
  thin = c(NA, setup$julia_thin),
  n_replications = c(setup$n_rep, setup$n_rep),
  n_models = c(length(setup$julia_models), length(setup$julia_models)),
  mechanisms = paste(setup$mechanisms, collapse = ";"),
  stringsAsFactors = FALSE
)
safe_write_csv(runtime_context, file.path(rev_root, "07_runtime", "runtime_context.csv"))




# 8. Optional summaries using existing helper functions------------------------

# These calls are kept optional because the exact column names may differ across
# intermediate versions of the code. The central raw files above are the source
# of truth and can be summarized in a separate table-generation script.

if (exists("harmonise_estimates_for_summary") && exists("summarise_estimates_over_reps")) {
  est_long <- harmonise_estimates_for_summary(
    ml_estimates = res_ml$ml_estimates,
    bayes_estimates = bayes_estimates
  )
  est_summary <- summarise_estimates_over_reps(est_long)
  safe_write_csv(est_long, file.path(rev_root, "01_parameter_recovery", "estimates_by_replication_M2.csv"))
  safe_write_csv(est_summary, file.path(rev_root, "01_parameter_recovery", "param_recovery_M2_all_mechanisms.csv"))
}

if (exists("make_bayes_factor_table") && !is.null(bayes_evidence)) {
  ev_mean <- make_bayes_factor_table(
    bayes_evidence,
    reference_model = "M2_true",
    fixed_point = "mean"
  )
  safe_write_csv(ev_mean, file.path(rev_root, "03_evidence", "evidence_vs_M2_all_mechanisms.csv"))
}

if (exists("summarise_predictive_criteria") && !is.null(bayes_pred)) {
  pred_summary <- summarise_predictive_criteria(bayes_pred)
  safe_write_csv(pred_summary, file.path(rev_root, "04_predictive_criteria", "predictive_criteria_all_workflows.csv"))
}

message("Core workflow finished. Results written to: ", rev_root)


# 9. Post-processing ----------------------------------------------------------

source("R/12_harmonise_outputs.R")
source("R/13_export_metadata.R")
source("R/18_postprocess_diagnostics.R")
source("R/14_validate_outputs.R")
source("R/19_make_paper_tables_figures.R")

harmonise_outputs(rev_out)

export_missingness_metadata(
  index_file = file.path(data_root, "simulation_index.csv"),
  out_dir = file.path(rev_out, "00_metadata")
)

export_complete_case_rates(
  index_file = file.path(data_root, "simulation_index.csv"),
  out_dir = file.path(rev_out, "00_metadata")
)

postprocess_reviewer_diagnostics(
  root = rev_out,
  make_plots = TRUE,
  trace_terms = c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a", "loglik")
)

validate_outputs(rev_out)




# 10. Plotting, tables, finale results ----------------------------------------

# Table helpers already loaded above

make_paper_tables_figures(
  root = rev_out,
  trace_terms = c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a", "loglik"),
  make_traceplots_flag = TRUE
)

readr::read_csv(
  file.path(rev_out, "08_paper_tables", "table_model_selection_rates_compact.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(criterion %in% c("BIC", "log_marginal_likelihood", "WAIC")) |>
  dplyr::arrange(mechanism, criterion, workflow, dplyr::desc(selection_rate)) |>
  print(n=Inf)


readr::read_csv(
  file.path(rev_out, "08_paper_tables", "table_runtime_compact.csv"),
  show_col_types = FALSE
) |>
  dplyr::distinct(source, workflow, mechanism)




###############
# Table helpers already loaded above

sel <- readr::read_csv(
  file.path(rev_out, "08_paper_tables", "table_model_selection_rates_compact.csv"),
  show_col_types = FALSE
)

sel_fixed <- sel |>
  dplyr::group_by(workflow, criterion, mechanism) |>
  dplyr::mutate(
    n_total = sum(n_selected, na.rm = TRUE),
    selection_rate = n_selected / n_total
  ) |>
  dplyr::ungroup()

readr::write_csv(
  sel_fixed,
  file.path(rev_out, "08_paper_tables", "table_model_selection_rates_compact.csv")
)


make_paper_tables_figures(
  root = rev_out,
  trace_terms = c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a", "loglik"),
  make_traceplots_flag = TRUE
)

readr::read_csv(
  file.path(rev_out, "08_paper_tables", "table_mcmc_diagnostics_compact.csv"),
  show_col_types = FALSE
)



# fix some 
par <- readr::read_csv(
  file.path(rev_out, "08_paper_tables", "table_param_recovery_M2_compact.csv"),
  show_col_types = FALSE
)

par_fixed <- par |>
  dplyr::arrange(workflow, mechanism, term, dplyr::desc(!is.na(coverage))) |>
  dplyr::group_by(workflow, mechanism, term) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

readr::write_csv(
  par_fixed,
  file.path(rev_out, "08_paper_tables", "table_param_recovery_M2_compact.csv")
)

message("Full workflow including post-processing finished: ", rev_root)
}
