
# Harmonise Bayesian and ML outputs after completed simulation runs.
# This script does not rerun any estimation. It only reads existing CSV files in
# FULL_OUT and creates derived files used by validation and paper-table scripts.

harmonise_workflow <- function(x) {
  dplyr::case_when(
    x %in% c("CART-BB medoid + conditional Bayesian MCMC (Julia)",
             "CART-BB medoid + conditional Bayesian MCMC") ~
      "CART-BB medoid + conditional Bayesian MCMC",
    x %in% c("Bayesian full-data MCMC (Julia)",
             "Bayesian full-data MCMC",
             "Before deletion + Bayesian MCMC") ~
      "Before deletion + Bayesian MCMC",
    x %in% c("Complete case + Bayesian MCMC (Julia)",
             "Complete case + Bayesian MCMC") ~
      "Complete case + Bayesian MCMC",
    x %in% c("Bayesian CART-BB MCMC (Julia)",
             "Bayesian CART-assisted MCMC (Julia)",
             "Bayesian CART-BB MCMC") ~
      "Bayesian CART-BB MCMC",
    x %in% c("MICE + ML", "MICE-CART + ML") ~
      "MICE-CART + ML",
    TRUE ~ as.character(x)
  )
}

canonicalise_parameter_term <- function(x) {
  out <- trimws(as.character(x))
  key <- tolower(gsub("[()_[:space:]]", "", out))
  out[key %in% c("intercept", "const", "cons")] <- "(Intercept)"
  out
}

ensure_workflow <- function(df) {
  if (!"workflow" %in% names(df)) {
    df$workflow <- NA_character_
  }
  if ("method" %in% names(df)) {
    miss <- is.na(df$workflow) | df$workflow == ""
    df$workflow[miss] <- df$method[miss]
  }
  df$workflow <- harmonise_workflow(df$workflow)
  df
}

default_true_values <- function() {
  c(
    "(Intercept)" = 0.20,
    "x2" = -1.50,
    "x3" = 0.80,
    "x4" = -0.30,
    "x5B" = -0.40,
    "x5C" = -0.30,
    "x7" = 0.00,
    "x8B" = 0.00,
    "x8C" = 0.00,
    "x10" = 0.00,
    "sigma2_a" = 1.00
  )
}

.safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.require_columns <- function(df, cols, name) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0L) {
    stop(name, " is missing required columns: ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

.first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else x[[1L]]
}

.mcse_prop <- function(p, n) sqrt(p * (1 - p) / pmax(n, 1))

.make_estimates_by_replication <- function(root) {
  bayes_path <- file.path(root, "01_parameter_recovery", "bayes_estimates_all_workflows.csv")
  ml_path    <- file.path(root, "01_parameter_recovery", "ml_estimates_all_workflows.csv")

  bayes <- .safe_read_csv(bayes_path)
  ml    <- .safe_read_csv(ml_path)

  pieces <- list()

  if (!is.null(bayes) && nrow(bayes) > 0L) {
    .require_columns(
      bayes,
      c("method", "model", "rep", "mechanism", "term", "posterior_mean",
        "posterior_sd", "ci_lower", "ci_upper", "component"),
      basename(bayes_path)
    )
    bayes <- ensure_workflow(bayes)
    bayes$term <- canonicalise_parameter_term(bayes$term)
    pieces[["bayes"]] <- data.frame(
      estimator = "Bayesian",
      method = bayes$method,
      workflow = bayes$workflow,
      model = bayes$model,
      rep = bayes$rep,
      mechanism = bayes$mechanism,
      term = bayes$term,
      estimate = bayes$posterior_mean,
      std_error = bayes$posterior_sd,
      ci_lower = bayes$ci_lower,
      ci_upper = bayes$ci_upper,
      component = bayes$component,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(ml) && nrow(ml) > 0L) {
    .require_columns(
      ml,
      c("method", "model", "rep", "mechanism", "term", "estimate", "std_error", "component"),
      basename(ml_path)
    )
    ml <- ensure_workflow(ml)
    ml$term <- canonicalise_parameter_term(ml$term)
    ml$ci_lower <- ifelse(is.finite(ml$std_error), ml$estimate - 1.96 * ml$std_error, NA_real_)
    ml$ci_upper <- ifelse(is.finite(ml$std_error), ml$estimate + 1.96 * ml$std_error, NA_real_)
    pieces[["ml"]] <- data.frame(
      estimator = "ML",
      method = ml$method,
      workflow = ml$workflow,
      model = ml$model,
      rep = ml$rep,
      mechanism = ml$mechanism,
      term = ml$term,
      estimate = ml$estimate,
      std_error = ml$std_error,
      ci_lower = ml$ci_lower,
      ci_upper = ml$ci_upper,
      component = ml$component,
      stringsAsFactors = FALSE
    )
  }

  if (length(pieces) == 0L) {
    stop("No Bayesian or ML estimate files found in ", file.path(root, "01_parameter_recovery"))
  }

  dplyr::bind_rows(pieces)
}

.summarise_parameter_recovery <- function(est, true_values) {
  truth <- data.frame(
    term = names(true_values),
    true_value = as.numeric(true_values),
    stringsAsFactors = FALSE
  )

  est2 <- dplyr::left_join(est, truth, by = "term")
  est2 <- dplyr::filter(est2, !is.na(.data$true_value), is.finite(.data$estimate))
  est2 <- dplyr::mutate(
    est2,
    bias = .data$estimate - .data$true_value,
    abs_bias = abs(.data$bias),
    sq_error = .data$bias^2,
    covered = dplyr::if_else(
      is.finite(.data$ci_lower) & is.finite(.data$ci_upper),
      .data$ci_lower <= .data$true_value & .data$ci_upper >= .data$true_value,
      NA
    )
  )

  dplyr::summarise(
    dplyr::group_by(est2, .data$estimator, .data$method, .data$workflow,
                    .data$mechanism, .data$model, .data$term, .data$component, .data$true_value),
    n_rep = dplyr::n_distinct(.data$rep),
    mean_estimate = mean(.data$estimate, na.rm = TRUE),
    sd_estimate = stats::sd(.data$estimate, na.rm = TRUE),
    bias = mean(.data$bias, na.rm = TRUE),
    mcse_bias = stats::sd(.data$bias, na.rm = TRUE) / sqrt(dplyr::n()),
    mean_abs_bias = mean(.data$abs_bias, na.rm = TRUE),
    rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)),
    coverage = mean(.data$covered, na.rm = TRUE),
    .groups = "drop"
  )
}

.make_selected_models <- function(root) {
  bayes_path <- file.path(root, "03_evidence", "evidence_by_replication.csv")
  ml_path    <- file.path(root, "02_model_selection", "ml_model_results_all_workflows.csv")

  bayes <- .safe_read_csv(bayes_path)
  ml    <- .safe_read_csv(ml_path)

  pieces <- list()

  if (!is.null(bayes) && nrow(bayes) > 0L) {
    .require_columns(
      bayes,
      c("method", "workflow", "model", "rep", "mechanism", "fixed_point",
        "log_marginal_likelihood"),
      basename(bayes_path)
    )
    bayes <- ensure_workflow(bayes)
    bayes <- dplyr::filter(bayes, .data$fixed_point == "mean",
                           is.finite(.data$log_marginal_likelihood))
    pieces[["bayes"]] <- bayes |>
      dplyr::group_by(.data$method, .data$workflow, .data$mechanism, .data$rep) |>
      dplyr::slice_max(order_by = .data$log_marginal_likelihood, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        estimator = "Bayesian",
        method = .data$method,
        workflow = .data$workflow,
        rep = .data$rep,
        mechanism = .data$mechanism,
        criterion = "log_marginal_likelihood",
        fixed_point = .data$fixed_point,
        selected_model = .data$model,
        selected_value = .data$log_marginal_likelihood
      )
  }

  if (!is.null(ml) && nrow(ml) > 0L) {
    .require_columns(
      ml,
      c("method", "model", "rep", "mechanism", "BIC"),
      basename(ml_path)
    )
    ml <- ensure_workflow(ml)
    ml <- dplyr::filter(ml, is.finite(.data$BIC))
    pieces[["ml"]] <- ml |>
      dplyr::group_by(.data$method, .data$workflow, .data$mechanism, .data$rep) |>
      dplyr::slice_min(order_by = .data$BIC, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        estimator = "ML",
        method = .data$method,
        workflow = .data$workflow,
        rep = .data$rep,
        mechanism = .data$mechanism,
        criterion = "BIC",
        fixed_point = NA_character_,
        selected_model = .data$model,
        selected_value = .data$BIC
      )
  }

  if (length(pieces) == 0L) {
    stop("No Bayesian evidence or ML model-result files found for model selection.")
  }

  dplyr::bind_rows(pieces)
}

.summarise_selection_rates <- function(sel) {
  sel |>
    dplyr::group_by(.data$estimator, .data$method, .data$workflow,
                    .data$mechanism, .data$criterion, .data$selected_model) |>
    dplyr::summarise(n_selected = dplyr::n(), .groups = "drop_last") |>
    dplyr::mutate(
      n_rep = sum(.data$n_selected),
      selection_rate = .data$n_selected / .data$n_rep,
      mcse_selection_rate = .mcse_prop(.data$selection_rate, .data$n_rep)
    ) |>
    dplyr::ungroup()
}

.make_runtime_tables <- function(root) {
  bayes_path <- file.path(root, "02_model_selection", "bayes_model_results_all_workflows.csv")
  ml_path    <- file.path(root, "02_model_selection", "ml_model_results_all_workflows.csv")
  comp_path  <- file.path(root, "07_runtime", "ml_runtime_components.csv")

  bayes <- .safe_read_csv(bayes_path)
  ml    <- .safe_read_csv(ml_path)
  comp  <- .safe_read_csv(comp_path)

  model_runtime_pieces <- list()

  if (!is.null(bayes) && nrow(bayes) > 0L && "runtime_seconds" %in% names(bayes)) {
    bayes <- ensure_workflow(bayes)
    model_runtime_pieces[["bayes"]] <- bayes |>
      dplyr::transmute(
        estimator = "Bayesian",
        method = .data$method,
        workflow = .data$workflow,
        mechanism = .data$mechanism,
        rep = .data$rep,
        model = .data$model,
        runtime_seconds = .data$runtime_seconds
      )
  }

  if (!is.null(ml) && nrow(ml) > 0L && "runtime_seconds" %in% names(ml)) {
    ml <- ensure_workflow(ml)
    model_runtime_pieces[["ml"]] <- ml |>
      dplyr::transmute(
        estimator = "ML",
        method = .data$method,
        workflow = .data$workflow,
        mechanism = .data$mechanism,
        rep = .data$rep,
        model = .data$model,
        runtime_seconds = .data$runtime_seconds
      )
  }

  runtime_model_long <- dplyr::bind_rows(model_runtime_pieces)

  runtime_by_model <- runtime_model_long |>
    dplyr::group_by(.data$estimator, .data$method, .data$workflow,
                    .data$mechanism, .data$model) |>
    dplyr::summarise(
      n_rep = dplyr::n_distinct(.data$rep),
      mean_runtime_seconds = mean(.data$runtime_seconds, na.rm = TRUE),
      median_runtime_seconds = stats::median(.data$runtime_seconds, na.rm = TRUE),
      sd_runtime_seconds = stats::sd(.data$runtime_seconds, na.rm = TRUE),
      min_runtime_seconds = min(.data$runtime_seconds, na.rm = TRUE),
      max_runtime_seconds = max(.data$runtime_seconds, na.rm = TRUE),
      .groups = "drop"
    )

  # Workflow-level runtimes. For ML imputation workflows, use component timings where available
  # to avoid counting the same common imputation once per candidate model.
  workflow_pieces <- list()

  if (!is.null(bayes) && nrow(bayes) > 0L && "runtime_seconds" %in% names(bayes)) {
    workflow_pieces[["bayes"]] <- bayes |>
      dplyr::group_by(.data$method, .data$workflow, .data$mechanism, .data$rep) |>
      dplyr::summarise(
        estimator = "Bayesian",
        total_runtime_seconds = sum(.data$runtime_seconds, na.rm = TRUE),
        .groups = "drop"
      )
  }

  if (!is.null(comp) && nrow(comp) > 0L) {
    .require_columns(
      comp,
      c("method", "rep", "mechanism", "model", "imputation_seconds", "model_fit_seconds"),
      basename(comp_path)
    )
    comp <- ensure_workflow(comp)
    comp_workflow <- comp |>
      dplyr::group_by(.data$method, .data$workflow, .data$mechanism, .data$rep) |>
      dplyr::summarise(
        estimator = "ML",
        total_runtime_seconds =
          max(.data$imputation_seconds, na.rm = TRUE) +
          sum(.data$model_fit_seconds, na.rm = TRUE),
        .groups = "drop"
      )
    workflow_pieces[["ml_components"]] <- comp_workflow
  }

  if (!is.null(ml) && nrow(ml) > 0L && "runtime_seconds" %in% names(ml)) {
    ml <- ensure_workflow(ml)
    ml_other <- ml |>
      dplyr::filter(!.data$method %in% c("MICE-CART + ML", "missForest + ML")) |>
      dplyr::group_by(.data$method, .data$workflow, .data$mechanism, .data$rep) |>
      dplyr::summarise(
        estimator = "ML",
        total_runtime_seconds = sum(.data$runtime_seconds, na.rm = TRUE),
        .groups = "drop"
      )
    workflow_pieces[["ml_other"]] <- ml_other
  }

  workflow_long <- dplyr::bind_rows(workflow_pieces)

  runtime_by_workflow_mechanism <- workflow_long |>
    dplyr::group_by(.data$estimator, .data$method, .data$workflow, .data$mechanism) |>
    dplyr::summarise(
      n_rep = dplyr::n_distinct(.data$rep),
      mean_total_runtime_seconds = mean(.data$total_runtime_seconds, na.rm = TRUE),
      median_total_runtime_seconds = stats::median(.data$total_runtime_seconds, na.rm = TRUE),
      sd_total_runtime_seconds = stats::sd(.data$total_runtime_seconds, na.rm = TRUE),
      min_total_runtime_seconds = min(.data$total_runtime_seconds, na.rm = TRUE),
      max_total_runtime_seconds = max(.data$total_runtime_seconds, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    runtime_by_model = runtime_by_model,
    runtime_by_workflow_mechanism = runtime_by_workflow_mechanism
  )
}

harmonise_outputs <- function(root,
                                       true_values = default_true_values(),
                                       overwrite = TRUE) {
  if (!dir.exists(root)) stop("Output root does not exist: ", root)

  dir.create(file.path(root, "01_parameter_recovery"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "02_model_selection"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "07_runtime"), recursive = TRUE, showWarnings = FALSE)

  estimates <- .make_estimates_by_replication(root)
  estimates_m2 <- dplyr::filter(estimates, .data$model == "M2_true")

  utils::write.csv(
    estimates,
    file.path(root, "01_parameter_recovery", "estimates_by_replication_all_models.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    estimates_m2,
    file.path(root, "01_parameter_recovery", "estimates_by_replication_M2.csv"),
    row.names = FALSE
  )

  param_all <- .summarise_parameter_recovery(estimates, true_values)
  param_m2 <- dplyr::filter(param_all, .data$model == "M2_true")

  utils::write.csv(
    param_all,
    file.path(root, "01_parameter_recovery", "param_recovery_all_models.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    param_m2,
    file.path(root, "01_parameter_recovery", "param_recovery_M2_all_mechanisms.csv"),
    row.names = FALSE
  )

  selected <- .make_selected_models(root)
  selection_rates <- .summarise_selection_rates(selected)

  utils::write.csv(
    selected,
    file.path(root, "02_model_selection", "selected_models_by_replication.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    selection_rates,
    file.path(root, "02_model_selection", "model_selection_rates_all_mechanisms.csv"),
    row.names = FALSE
  )

  rt <- .make_runtime_tables(root)
  utils::write.csv(
    rt$runtime_by_workflow_mechanism,
    file.path(root, "07_runtime", "runtime_by_workflow_mechanism.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    rt$runtime_by_model,
    file.path(root, "07_runtime", "runtime_by_model.csv"),
    row.names = FALSE
  )

  message("Harmonised outputs written to: ", normalizePath(root, winslash = "/", mustWork = FALSE))

  invisible(list(
    estimates_by_replication_all_models = estimates,
    estimates_by_replication_M2 = estimates_m2,
    param_recovery_all_models = param_all,
    param_recovery_M2_all_mechanisms = param_m2,
    selected_models_by_replication = selected,
    model_selection_rates_all_mechanisms = selection_rates,
    runtime_by_workflow_mechanism = rt$runtime_by_workflow_mechanism,
    runtime_by_model = rt$runtime_by_model
  ))
}
