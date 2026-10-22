# -----------------------------------------------------------------------------
# 17_sensitivity_summaries.R
# Fixed-point and prior-sensitivity summaries for Bayesian evidence measures.
#
# Reviewer-relevant outputs:
#   05_sensitivity/fixed_point_sensitivity.csv
#   05_sensitivity/fixed_point_selected_models.csv
#   05_sensitivity/prior_sensitivity_tau2.csv       (if prior_tau2 is available)
#   05_sensitivity/prior_selected_models_tau2.csv   (if prior_tau2 is available)
# -----------------------------------------------------------------------------

harmonise_workflow_label_for_sensitivity <- function(x) {
  dplyr::case_when(
    x %in% c("Bayesian full-data MCMC (Julia)", "Bayesian full-data MCMC") ~
      "Before deletion + Bayesian MCMC",
    x %in% c("Complete case + Bayesian MCMC (Julia)", "Complete case + Bayesian MCMC") ~
      "Complete case + Bayesian MCMC",
    x %in% c("Bayesian CART-BB MCMC (Julia)", "Bayesian CART-assisted MCMC (Julia)",
             "Bayesian CART-BB MCMC") ~
      "Bayesian CART-BB MCMC",
    TRUE ~ x
  )
}

make_fixed_point_sensitivity <- function(root,
                                         evidence_file = file.path(root, "03_evidence", "evidence_by_replication.csv")) {
  out_dir <- file.path(root, "05_sensitivity")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(evidence_file)) {
    warning("Evidence file not found: ", evidence_file)
    return(invisible(NULL))
  }

  ev <- readr::read_csv(evidence_file, show_col_types = FALSE)

  if (!"workflow" %in% names(ev)) ev$workflow <- ev$method
  if (!"fixed_point" %in% names(ev)) {
    stop("Evidence file must contain a 'fixed_point' column with values such as mean, median, high_density.")
  }

  ev <- ev |>
    dplyr::mutate(
      workflow = harmonise_workflow_label_for_sensitivity(.data$workflow),
      fixed_point = as.character(.data$fixed_point),
      mechanism = toupper(.data$mechanism)
    )

  fp <- ev |>
    dplyr::group_by(.data$workflow, .data$mechanism, .data$fixed_point, .data$model) |>
    dplyr::summarise(
      mean_log_evidence = mean(.data$log_marginal_likelihood, na.rm = TRUE),
      sd_log_evidence = stats::sd(.data$log_marginal_likelihood, na.rm = TRUE),
      median_log_evidence = stats::median(.data$log_marginal_likelihood, na.rm = TRUE),
      n_rep = dplyr::n_distinct(.data$rep),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$workflow, .data$mechanism, .data$fixed_point) |>
    dplyr::mutate(
      rank = rank(-.data$mean_log_evidence, ties.method = "min"),
      selected_model = .data$model[which.max(.data$mean_log_evidence)]
    ) |>
    dplyr::ungroup()

  fp_selected <- fp |>
    dplyr::filter(.data$rank == 1L) |>
    dplyr::select(.data$workflow, .data$mechanism, .data$fixed_point, .data$selected_model,
                  .data$mean_log_evidence, .data$sd_log_evidence, .data$n_rep)

  readr::write_csv(fp, file.path(out_dir, "fixed_point_sensitivity.csv"))
  readr::write_csv(fp_selected, file.path(out_dir, "fixed_point_selected_models.csv"))

  message("Fixed-point sensitivity written to: ", out_dir)
  invisible(list(fixed_point_sensitivity = fp, selected = fp_selected))
}

make_prior_sensitivity_tau2 <- function(root,
                                        prior_file = file.path(root, "05_sensitivity", "prior_evidence_by_replication.csv")) {
  out_dir <- file.path(root, "05_sensitivity")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(prior_file)) {
    warning("Prior-sensitivity input not found: ", prior_file,
            "\nExpected columns include: rep, mechanism, workflow/method, model, prior_tau2, log_marginal_likelihood.")
    return(invisible(NULL))
  }

  pr <- readr::read_csv(prior_file, show_col_types = FALSE)
  if (!"workflow" %in% names(pr)) pr$workflow <- pr$method
  if (!"prior_tau2" %in% names(pr)) {
    stop("prior_file must contain a 'prior_tau2' column.")
  }

  pr <- pr |>
    dplyr::mutate(
      workflow = harmonise_workflow_label_for_sensitivity(.data$workflow),
      mechanism = toupper(.data$mechanism),
      prior_tau2 = as.numeric(.data$prior_tau2)
    )

  sens <- pr |>
    dplyr::group_by(.data$workflow, .data$mechanism, .data$prior_tau2, .data$model) |>
    dplyr::summarise(
      mean_log_evidence = mean(.data$log_marginal_likelihood, na.rm = TRUE),
      sd_log_evidence = stats::sd(.data$log_marginal_likelihood, na.rm = TRUE),
      median_log_evidence = stats::median(.data$log_marginal_likelihood, na.rm = TRUE),
      n_rep = dplyr::n_distinct(.data$rep),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$workflow, .data$mechanism, .data$prior_tau2) |>
    dplyr::mutate(
      rank = rank(-.data$mean_log_evidence, ties.method = "min"),
      selected_model = .data$model[which.max(.data$mean_log_evidence)]
    ) |>
    dplyr::ungroup()

  selected <- sens |>
    dplyr::filter(.data$rank == 1L) |>
    dplyr::select(.data$workflow, .data$mechanism, .data$prior_tau2, .data$selected_model,
                  .data$mean_log_evidence, .data$sd_log_evidence, .data$n_rep)

  readr::write_csv(sens, file.path(out_dir, "prior_sensitivity_tau2.csv"))
  readr::write_csv(selected, file.path(out_dir, "prior_selected_models_tau2.csv"))

  message("Prior-sensitivity summaries written to: ", out_dir)
  invisible(list(prior_sensitivity = sens, selected = selected))
}
