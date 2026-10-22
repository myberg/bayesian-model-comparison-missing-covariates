# -----------------------------------------------------------------------------
# 19_make_paper_tables_figures.R
#
# Build compact reviewer-facing tables and figures from the CSV outputs.
#
# Usage from the analysis/ project root:
#   source("R/19_make_paper_tables_figures.R")
#   make_paper_tables_figures(root = "outputs/smoke10")
#
# For the final run, replace smoke10 by the final run folder, e.g.
#   make_paper_tables_figures(root = "outputs/full_20260616_0900")
# -----------------------------------------------------------------------------

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

safe_max_abs <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(abs(x))
}


safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}



# ---- Package handling --------------------------------------------------------
require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      "\nInstall them first, e.g.: install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), "))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_packages(c("readr", "dplyr", "tidyr", "stringr", "purrr", "ggplot2"))

# ---- General helpers ---------------------------------------------------------
message2 <- function(...) message(sprintf(...))

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

path_first_existing <- function(root, candidates) {
  paths <- file.path(root, candidates)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1L]]
}

safe_read_csv <- function(path, required = FALSE) {
  if (is.na(path) || !file.exists(path)) {
    if (required) stop("Required file not found: ", path, call. = FALSE)
    message2("Skipping missing file: %s", path)
    return(NULL)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0L) {
    if (required) stop("Required file is empty: ", path, call. = FALSE)
    message2("Skipping empty file: %s", path)
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

coalesce_col <- function(dat, candidates, default = NA) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0L) return(rep(default, nrow(dat)))
  dat[[hit[[1L]]]]
}

canonicalise_parameter_term <- function(x) {
  out <- trimws(as.character(x))
  key <- tolower(gsub("[()_[:space:]]", "", out))
  out[key %in% c("intercept", "const", "cons")] <- "(Intercept)"
  out
}

ensure_workflow <- function(dat) {
  if (is.null(dat)) return(NULL)
  if (!"workflow" %in% names(dat) && "method" %in% names(dat)) dat$workflow <- dat$method
  if ("workflow" %in% names(dat)) {
    wf <- as.character(dat$workflow)
    if ("method" %in% names(dat)) {
      method <- as.character(dat$method)
      is_conditional_medoid <- stringr::str_detect(
        method,
        stringr::regex("medoid.*conditional|conditional.*medoid", ignore_case = TRUE)
      )
      wf[is_conditional_medoid] <- method[is_conditional_medoid]
    }
    source_col <- if ("source" %in% names(dat)) as.character(dat$source) else rep(NA_character_, nrow(dat))
    
    dat$workflow <- dplyr::case_when(
      # Put exact/specific ML labels before generic before-deletion rules.
      stringr::str_detect(wf, stringr::regex("Before deletion \\+ ML", ignore_case = TRUE)) ~ "Before deletion + ML",
      stringr::str_detect(wf, stringr::regex("Complete case \\+ ML", ignore_case = TRUE)) ~ "Complete case + ML",
      stringr::str_detect(wf, stringr::regex("MICE", ignore_case = TRUE)) ~ "MICE-CART + ML",
      stringr::str_detect(wf, stringr::regex("missForest", ignore_case = TRUE)) ~ "missForest + ML",
      
      # Bayesian workflows. Some Julia outputs use full.data/full as the data-version label.
      stringr::str_detect(wf, stringr::regex("medoid.*conditional|conditional.*medoid", ignore_case = TRUE)) ~ "CART-BB medoid + conditional Bayesian MCMC",
      stringr::str_detect(wf, stringr::regex("CART.BB|CART-BB", ignore_case = TRUE)) ~ "Bayesian CART-BB MCMC",
      stringr::str_detect(wf, stringr::regex("Complete case.*Bayesian|Complete case \\+ Bayesian", ignore_case = TRUE)) ~ "Complete case + Bayesian MCMC",
      stringr::str_detect(wf, stringr::regex("Before deletion.*Bayesian|Before deletion \\+ Bayesian|full.data", ignore_case = TRUE)) ~ "Before deletion + Bayesian MCMC",
      stringr::str_detect(wf, stringr::regex("before deletion", ignore_case = TRUE)) &
        stringr::str_detect(source_col, stringr::regex("ML", ignore_case = TRUE)) ~ "Before deletion + ML",
      stringr::str_detect(wf, stringr::regex("before deletion", ignore_case = TRUE)) ~ "Before deletion + Bayesian MCMC",
      
      TRUE ~ wf
    )
  }
  dat
}

round_numeric <- function(dat, digits = 3) {
  if (is.null(dat)) return(NULL)
  dplyr::mutate(dat, dplyr::across(where(is.numeric), ~ round(.x, digits)))
}

escape_tex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

write_latex_table <- function(dat, path, caption = NULL, label = NULL, digits = 3) {
  dat <- as.data.frame(round_numeric(dat, digits = digits))
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  
  if (requireNamespace("knitr", quietly = TRUE)) {
    tex <- knitr::kable(
      dat,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      caption = caption,
      label = label,
      escape = TRUE,
      align = "l"
    )
    writeLines(tex, path, useBytes = TRUE)
    return(invisible(path))
  }
  
  # Minimal fallback if knitr is not installed.
  header <- paste(escape_tex(names(dat)), collapse = " & ")
  rows <- apply(dat, 1, function(r) paste(escape_tex(r), collapse = " & "))
  out <- c(
    if (!is.null(caption)) paste0("% ", caption) else NULL,
    if (!is.null(label)) paste0("% ", label) else NULL,
    "\\begin{tabular}{", paste(rep("l", ncol(dat)), collapse = ""), "}",
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  writeLines(out, path, useBytes = TRUE)
  invisible(path)
}

write_table_outputs <- function(dat, name, out_tables, caption = NULL, label = NULL, digits = 3) {
  if (is.null(dat) || nrow(dat) == 0L) {
    message2("Skipping empty table: %s", name)
    return(invisible(NULL))
  }
  make_dir(out_tables)
  csv_path <- file.path(out_tables, paste0(name, ".csv"))
  tex_path <- file.path(out_tables, paste0(name, ".tex"))
  readr::write_csv(round_numeric(dat, digits = digits), csv_path)
  write_latex_table(dat, tex_path, caption = caption, label = label, digits = digits)
  message2("Wrote table: %s", csv_path)
  invisible(list(csv = csv_path, tex = tex_path))
}

save_plot_outputs <- function(plot, name, out_figures, width = 8.0, height = 5.0) {
  if (is.null(plot)) return(invisible(NULL))
  make_dir(out_figures)
  png_path <- file.path(out_figures, paste0(name, ".png"))
  pdf_path <- file.path(out_figures, paste0(name, ".pdf"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 320)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height)
  message2("Wrote figure: %s", png_path)
  invisible(list(png = png_path, pdf = pdf_path))
}

base_theme <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

workflow_shape_values <- c(
  "Bayesian CART-BB MCMC" = 16,
  "Before deletion + Bayesian MCMC" = 17,
  "Complete case + Bayesian MCMC" = 15,
  "Before deletion + ML" = 4,
  "Complete case + ML" = 3,
  "MICE-CART + ML" = 7,
  "missForest + ML" = 8
)

workflow_shape_scale <- function() {
  ggplot2::scale_shape_manual(
    values = workflow_shape_values,
    breaks = names(workflow_shape_values),
    na.translate = FALSE
  )
}

# ---- Table builders ----------------------------------------------------------
make_missingness_table <- function(root, out_tables) {
  path <- path_first_existing(root, c("00_metadata/realized_missingness_rates.csv"))
  dat <- safe_read_csv(path)
  if (is.null(dat)) return(NULL)
  
  var_col <- intersect(c("variable", "term", "covariate"), names(dat))[1]
  if (is.na(var_col)) return(NULL)
  
  rate_col <- intersect(c("mean_missing_rate", "missing_rate.mean", "missing_rate", "rate"), names(dat))[1]
  if (is.na(rate_col)) {
    num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
    rate_col <- num_cols[1]
  }
  
  out <- dat |>
    dplyr::mutate(
      mechanism = as.character(.data$mechanism),
      variable = as.character(.data[[var_col]]),
      missing_rate = as.numeric(.data[[rate_col]])
    ) |>
    dplyr::group_by(.data$mechanism, .data$variable) |>
    dplyr::summarise(missing_rate = mean(.data$missing_rate, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data$mechanism, values_from = .data$missing_rate) |>
    dplyr::arrange(.data$variable)
  
  write_table_outputs(
    out, "table_missingness_rates_compact", out_tables,
    caption = "Realised missingness rates by mechanism and covariate.",
    label = "tab:missingness_rates"
  )
  out
}

make_complete_case_table <- function(root, out_tables) {
  path <- path_first_existing(root, c("00_metadata/complete_case_rates.csv"))
  dat <- safe_read_csv(path)
  if (is.null(dat)) return(NULL)
  
  rate_col <- intersect(c("mean_complete_case_rate", "complete_case_rate.mean", "complete_case_rate", "rate"), names(dat))[1]
  if (is.na(rate_col)) {
    num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
    rate_col <- num_cols[1]
  }
  
  out <- dat |>
    dplyr::mutate(
      mechanism = as.character(.data$mechanism),
      complete_case_rate = as.numeric(.data[[rate_col]])
    ) |>
    dplyr::group_by(.data$mechanism) |>
    dplyr::summarise(complete_case_rate = mean(.data$complete_case_rate, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$mechanism)
  
  write_table_outputs(
    out, "table_complete_case_rates", out_tables,
    caption = "Average complete-case rates by missing-data mechanism.",
    label = "tab:complete_case_rates"
  )
  out
}

make_param_recovery_table <- function(root, out_tables, terms_keep = NULL) {
  path <- path_first_existing(root, c(
    "01_parameter_recovery/param_recovery_M2_all_mechanisms.csv",
    "01_parameter_recovery/param_recovery_M2.csv"
  ))
  dat <- safe_read_csv(path)
  dat <- ensure_workflow(dat)
  if (is.null(dat)) return(NULL)
  
  if (is.null(terms_keep)) {
    terms_keep <- c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a")
  }
  
  out <- dat |>
    dplyr::mutate(
      term = canonicalise_parameter_term(.data$term),
      mean_estimate = as.numeric(coalesce_col(dplyr::cur_data_all(), c("mean_estimate", "estimate_mean", "posterior_mean"))),
      bias = as.numeric(coalesce_col(dplyr::cur_data_all(), c("bias", "mean_bias"))),
      rmse = as.numeric(coalesce_col(dplyr::cur_data_all(), c("rmse", "RMSE"))),
      mcse_bias = as.numeric(coalesce_col(dplyr::cur_data_all(), c("mcse_bias"))),
      mcse_rmse = as.numeric(coalesce_col(dplyr::cur_data_all(), c("mcse_rmse"))),
      coverage = as.numeric(coalesce_col(dplyr::cur_data_all(), c("coverage", "coverage_95"))),
      mcse_coverage = as.numeric(coalesce_col(dplyr::cur_data_all(), c("mcse_coverage"))),
      n_rep = as.numeric(coalesce_col(dplyr::cur_data_all(), c("n_rep", "n")))
    ) |>
    dplyr::filter(.data$term %in% terms_keep) |>
    dplyr::mutate(term = factor(.data$term, levels = terms_keep)) |>
    dplyr::select(dplyr::any_of(c(
      "workflow", "mechanism", "term", "true_value", "mean_estimate", "bias",
      "mcse_bias", "rmse", "mcse_rmse", "coverage", "mcse_coverage", "n_rep"
    ))) |>
    dplyr::arrange(.data$mechanism, .data$workflow, factor(.data$term, levels = terms_keep))

  if (!"(Intercept)" %in% out$term) {
    stop("Parameter-recovery table is missing the intercept after term harmonisation.", call. = FALSE)
  }
  
  write_table_outputs(
    out, "table_param_recovery_M2_compact", out_tables,
    caption = "Parameter recovery for the data-generating model.",
    label = "tab:param_recovery_m2"
  )
  out
}

make_model_selection_table <- function(root, out_tables) {
  path <- path_first_existing(root, c("02_model_selection/model_selection_rates_all_mechanisms.csv"))
  dat <- safe_read_csv(path)
  dat <- ensure_workflow(dat)
  if (is.null(dat)) return(NULL)
  
  # The raw selection-rate file may contain either `n_selected` or simply `n`.
  # In some intermediate versions, `n_total` was copied row-wise from `n`, which
  # makes every selection_rate equal to 1. Recompute totals within
  # workflow x criterion x mechanism to obtain actual selection frequencies.
  out <- dat |>
    dplyr::mutate(
      selected_model = as.character(coalesce_col(dplyr::cur_data_all(), c("selected_model", "model"))),
      n_selected = as.numeric(coalesce_col(dplyr::cur_data_all(), c("n_selected", "count", "n", "Freq"))),
      selection_rate_raw = as.numeric(coalesce_col(dplyr::cur_data_all(), c("selection_rate", "rate")))
    ) |>
    dplyr::select(dplyr::any_of(c("workflow", "criterion", "mechanism", "selected_model", "n_selected", "selection_rate_raw"))) |>
    dplyr::group_by(.data$workflow, .data$criterion, .data$mechanism) |>
    dplyr::mutate(
      n_total = sum(.data$n_selected, na.rm = TRUE),
      selection_rate = dplyr::if_else(
        is.finite(.data$n_total) & .data$n_total > 0,
        .data$n_selected / .data$n_total,
        .data$selection_rate_raw
      ),
      mcse_selection_rate = sqrt(.data$selection_rate * (1 - .data$selection_rate) / .data$n_total)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::any_of(c("workflow", "criterion", "mechanism", "selected_model", "selection_rate", "mcse_selection_rate", "n_selected", "n_total"))) |>
    dplyr::arrange(.data$mechanism, .data$workflow, .data$criterion, dplyr::desc(.data$selection_rate))
  
  write_table_outputs(
    out, "table_model_selection_rates_compact", out_tables,
    caption = "Model-selection frequencies by workflow, mechanism and criterion.",
    label = "tab:model_selection_rates"
  )
  out
}

make_evidence_vs_ref <- function(ev, reference_model = "M2_true") {
  ev <- ensure_workflow(ev)
  if (is.null(ev) || !"log_marginal_likelihood" %in% names(ev)) return(NULL)
  
  keys <- intersect(c("workflow", "method", "rep", "mechanism", "data_version", "impute_mode", "fixed_point"), names(ev))
  if (!"workflow" %in% keys && "method" %in% keys) {
    ev$workflow <- ev$method
    keys <- unique(c("workflow", setdiff(keys, "method")))
  }
  ref <- ev |>
    dplyr::filter(.data$model == reference_model) |>
    dplyr::select(dplyr::all_of(keys), log_marginal_likelihood_ref = .data$log_marginal_likelihood)
  
  out <- ev |>
    dplyr::left_join(ref, by = keys) |>
    dplyr::mutate(
      reference_model = reference_model,
      lnBF_model_vs_ref = .data$log_marginal_likelihood - .data$log_marginal_likelihood_ref
    )
  out
}

make_evidence_table <- function(root, out_tables) {
  path <- path_first_existing(root, c(
    "03_evidence/evidence_vs_M2_all_mechanisms.csv",
    "03_evidence/evidence_by_replication.csv"
  ))
  dat <- safe_read_csv(path)
  dat <- ensure_workflow(dat)
  if (is.null(dat)) return(NULL)
  
  if (!"lnBF_model_vs_ref" %in% names(dat)) {
    dat <- make_evidence_vs_ref(dat, reference_model = "M2_true")
  }
  if (is.null(dat)) return(NULL)
  
  out <- dat |>
    dplyr::filter(is.na(.data$fixed_point) | .data$fixed_point == "mean") |>
    dplyr::group_by(.data$workflow, .data$mechanism, .data$model, .data$reference_model) |>
    dplyr::summarise(
      mean_lnBF_vs_M2 = mean(.data$lnBF_model_vs_ref, na.rm = TRUE),
      sd_lnBF_vs_M2 = stats::sd(.data$lnBF_model_vs_ref, na.rm = TRUE),
      median_lnBF_vs_M2 = stats::median(.data$lnBF_model_vs_ref, na.rm = TRUE),
      n_rep = dplyr::n_distinct(.data$rep),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$mechanism, .data$workflow, dplyr::desc(.data$mean_lnBF_vs_M2))
  
  write_table_outputs(
    out, "table_evidence_vs_M2_compact", out_tables,
    caption = "Bayesian evidence relative to the data-generating model using the posterior-mean fixed point.",
    label = "tab:evidence_vs_m2"
  )
  out
}

make_fixed_point_table <- function(root, out_tables) {
  path <- path_first_existing(root, c("05_sensitivity/fixed_point_sensitivity.csv"))
  dat <- safe_read_csv(path)
  dat <- ensure_workflow(dat)
  if (is.null(dat)) return(NULL)
  
  keep <- intersect(c("workflow", "mechanism", "fixed_point", "selected_model", "mean_log_evidence", "sd_log_evidence", "n_rep"), names(dat))
  out <- dat |>
    dplyr::select(dplyr::all_of(keep)) |>
    dplyr::arrange(.data$mechanism, .data$workflow, .data$fixed_point)
  
  write_table_outputs(
    out, "table_fixed_point_sensitivity", out_tables,
    caption = "Sensitivity of evidence-based model choice to the Chib evaluation point.",
    label = "tab:fixed_point_sensitivity"
  )
  out
}

make_mcmc_diagnostics_table <- function(root, out_tables) {
  ess_path <- path_first_existing(root, c("06_mcmc_diagnostics/ess_summary.csv"))
  ac_path <- path_first_existing(root, c("06_mcmc_diagnostics/autocorr_summary.csv"))
  mcse_path <- path_first_existing(root, c("06_mcmc_diagnostics/mcse_summary.csv"))
  
  ess <- safe_read_csv(ess_path)
  ac <- safe_read_csv(ac_path)
  mcse <- safe_read_csv(mcse_path)
  
  ess <- ensure_workflow(ess)
  ac <- ensure_workflow(ac)
  mcse <- ensure_workflow(mcse)
  
  if (is.null(ess)) return(NULL)
  
  keys <- intersect(c("workflow", "mechanism", "model", "rep", "chain", "term"), names(ess))
  
  out <- ess
  if (!is.null(ac)) {
    join_keys <- intersect(keys, names(ac))
    out <- dplyr::left_join(out, ac, by = join_keys, suffix = c("", "_ac"))
  }
  if (!is.null(mcse)) {
    join_keys <- intersect(keys, names(mcse))
    out <- dplyr::left_join(out, mcse, by = join_keys, suffix = c("", "_mcse"))
  }
  
  first_numeric_col <- function(dat, candidates) {
    hit <- candidates[candidates %in% names(dat)]
    if (length(hit) == 0L) return(rep(NA_real_, nrow(dat)))
    as.numeric(dat[[hit[[1L]]]])
  }
  
  # Different diagnostic helper versions used different column names.
  # The current diagnostic file uses ess_approx; older variants used
  # ess_bulk_approx or ess_bulk.
  out <- out |>
    dplyr::mutate(
      ess_bulk_for_summary = first_numeric_col(
        dplyr::cur_data_all(),
        c("ess_bulk", "ess_bulk_approx", "ess_approx")
      ),
      ess_tail_for_summary = first_numeric_col(
        dplyr::cur_data_all(),
        c("ess_tail", "ess_tail_approx")
      ),
      autocorr_lag1_for_summary = first_numeric_col(
        dplyr::cur_data_all(),
        c("autocorr_lag1", "autocorr_lag1_ac")
      ),
      mcse_for_summary = first_numeric_col(
        dplyr::cur_data_all(),
        c("mcse", "mcse_mean", "mcse_mean_mcse")
      )
    )
  
  out <- out |>
    dplyr::group_by(.data$workflow, .data$mechanism, .data$term) |>
    dplyr::summarise(
      min_ess_bulk = safe_min(.data$ess_bulk_for_summary),
      mean_ess_bulk = safe_mean(.data$ess_bulk_for_summary),
      mean_autocorr_lag1 = safe_mean(.data$autocorr_lag1_for_summary),
      mean_mcse = safe_mean(.data$mcse_for_summary),
      n_series = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$mechanism, .data$workflow, .data$term)
  
  write_table_outputs(
    out, "table_mcmc_diagnostics_compact", out_tables,
    caption = "MCMC diagnostic summaries for the diagnostic runs.",
    label = "tab:mcmc_diagnostics"
  )
  out
}

make_runtime_table <- function(root, out_tables) {
  path <- path_first_existing(root, c("07_runtime/runtime_by_workflow_mechanism.csv"))
  dat <- safe_read_csv(path)
  dat <- ensure_workflow(dat)
  if (is.null(dat)) return(NULL)

  # The harmonised runtime table may come from either the older summary
  # schema (source/runtime_seconds.*) or the revised full100 schema
  # (estimator/mean_total_runtime_seconds).  Keep whichever columns are
  # available and arrange only by columns that actually exist.
  candidate_cols <- c(
    "source", "estimator", "method", "workflow", "mechanism",
    "data_version", "impute_mode", "n_rep",
    "runtime_seconds.mean", "runtime_seconds.median",
    "runtime_seconds.q25", "runtime_seconds.q75", "runtime_seconds.sd",
    "mean_total_runtime_seconds", "median_total_runtime_seconds",
    "sd_total_runtime_seconds", "min_total_runtime_seconds",
    "max_total_runtime_seconds"
  )
  cols <- intersect(candidate_cols, names(dat))

  out <- dat |>
    dplyr::select(dplyr::all_of(cols))

  arrange_cols <- intersect(c("source", "estimator", "workflow", "mechanism"), names(out))
  if (length(arrange_cols) > 0L) {
    out <- out |> dplyr::arrange(dplyr::across(dplyr::all_of(arrange_cols)))
  }

  write_table_outputs(
    out, "table_runtime_compact", out_tables,
    caption = "Runtime by workflow and missing-data mechanism.",
    label = "tab:runtime"
  )
  out
}

# ---- Figure builders ---------------------------------------------------------
plot_model_selection_rates <- function(sel_rates) {
  if (is.null(sel_rates) || nrow(sel_rates) == 0L) return(NULL)
  if (!all(c("workflow", "criterion", "mechanism", "selected_model", "selection_rate") %in% names(sel_rates))) return(NULL)
  
  sel_rates |>
    dplyr::filter(!is.na(.data$selection_rate)) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$selected_model, y = .data$selection_rate, fill = .data$workflow)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(.data$mechanism ~ .data$criterion) +
    ggplot2::labs(x = "Selected model", y = "Selection rate", fill = "Workflow") +
    base_theme()
}

plot_parameter_bias <- function(param_tab) {
  if (is.null(param_tab) || nrow(param_tab) == 0L || !"bias" %in% names(param_tab)) return(NULL)
  
  param_tab |>
    ggplot2::ggplot(ggplot2::aes(x = .data$term, y = .data$bias, group = .data$workflow)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$workflow), position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::facet_wrap(~ mechanism, nrow = 1) +
    ggplot2::labs(x = "Parameter", y = "Bias", shape = "Workflow") +
    workflow_shape_scale() +
    base_theme()
}

plot_parameter_rmse <- function(param_tab) {
  if (is.null(param_tab) || nrow(param_tab) == 0L || !"rmse" %in% names(param_tab)) return(NULL)
  
  param_tab |>
    ggplot2::ggplot(ggplot2::aes(x = .data$term, y = .data$rmse, group = .data$workflow)) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$workflow), position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::facet_wrap(~ mechanism, nrow = 1) +
    ggplot2::labs(x = "Parameter", y = "RMSE", shape = "Workflow") +
    workflow_shape_scale() +
    base_theme()
}

plot_evidence_heatmap <- function(evidence_tab) {
  if (is.null(evidence_tab) || nrow(evidence_tab) == 0L || !"mean_lnBF_vs_M2" %in% names(evidence_tab)) return(NULL)
  
  evidence_tab |>
    ggplot2::ggplot(ggplot2::aes(x = .data$model, y = .data$workflow, fill = .data$mean_lnBF_vs_M2)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ mechanism) +
    ggplot2::labs(x = "Model", y = "Workflow", fill = "Mean lnBF vs M2") +
    base_theme()
}

plot_fixed_point_sensitivity <- function(fp_tab) {
  if (is.null(fp_tab) || nrow(fp_tab) == 0L) return(NULL)
  if (!all(c("workflow", "mechanism", "fixed_point", "mean_log_evidence") %in% names(fp_tab))) return(NULL)
  
  fp_tab |>
    ggplot2::ggplot(ggplot2::aes(x = .data$fixed_point, y = .data$mean_log_evidence, group = .data$workflow)) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$workflow), position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::facet_wrap(~ mechanism, scales = "free_y") +
    ggplot2::labs(x = "Evaluation point", y = "Mean log evidence", shape = "Workflow") +
    workflow_shape_scale() +
    base_theme()
}

plot_runtime <- function(runtime_tab) {
  if (is.null(runtime_tab) || nrow(runtime_tab) == 0L) return(NULL)
  y_col <- intersect(c("runtime_seconds.mean", "mean_runtime_seconds", "runtime_seconds"), names(runtime_tab))[1]
  if (is.na(y_col)) return(NULL)
  
  runtime_tab |>
    dplyr::mutate(runtime_minutes = as.numeric(.data[[y_col]]) / 60) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$workflow, y = .data$runtime_minutes, fill = .data$mechanism)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(x = "Workflow", y = "Mean runtime per model fit (minutes)", fill = "Mechanism") +
    base_theme()
}

plot_mcmc_ess <- function(mcmc_tab) {
  if (is.null(mcmc_tab) || nrow(mcmc_tab) == 0L || !"mean_ess_bulk" %in% names(mcmc_tab)) return(NULL)
  
  mcmc_tab |>
    dplyr::filter(is.finite(.data$mean_ess_bulk)) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$term, y = .data$mean_ess_bulk, shape = .data$workflow)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::facet_wrap(~ mechanism, scales = "free_y") +
    ggplot2::labs(x = "Parameter", y = "Mean approximate bulk ESS", shape = "Workflow") +
    workflow_shape_scale() +
    base_theme()
}

make_traceplots <- function(root, out_figures,
                            trace_terms = c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a", "loglik"),
                            max_files = Inf) {
  draw_dir <- file.path(root, "06_mcmc_diagnostics")
  files <- list.files(draw_dir, pattern = "^draws_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0L) {
    message2("No draw files found for traceplots in %s", draw_dir)
    return(invisible(NULL))
  }
  files <- head(files, max_files)
  
  draws <- purrr::map_dfr(files, ~ readr::read_csv(.x, show_col_types = FALSE, progress = FALSE))
  draws <- ensure_workflow(draws)
  draws <- draws |>
    dplyr::mutate(term = canonicalise_parameter_term(.data$term)) |>
    dplyr::filter(.data$term %in% trace_terms)
  
  if (nrow(draws) == 0L) {
    message2("No requested trace terms found in draw files.")
    return(invisible(NULL))
  }
  
  make_dir(file.path(out_figures, "traceplots"))
  
  # One compact plot per mechanism keeps the PDF readable.
  for (mech in sort(unique(draws$mechanism))) {
    dat <- draws |>
      dplyr::filter(.data$mechanism == mech)
    
    p <- dat |>
      ggplot2::ggplot(ggplot2::aes(x = .data$iter, y = .data$value)) +
      ggplot2::geom_line(linewidth = 0.25) +
      ggplot2::facet_grid(term ~ workflow, scales = "free_y") +
      ggplot2::labs(x = "Retained draw", y = "Value", title = paste("Traceplots:", mech)) +
      base_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5), legend.position = "none")
    
    save_plot_outputs(
      p,
      paste0("traceplots/fig_traceplots_", stringr::str_to_lower(mech)),
      out_figures,
      width = 10,
      height = max(5, 1.1 * length(unique(dat$term)))
    )
  }
  
  invisible(draws)
}

# ---- Main driver -------------------------------------------------------------
make_paper_tables_figures <- function(root = "outputs/smoke10",
                                      out_tables = file.path(root, "08_paper_tables"),
                                      out_figures = file.path(root, "09_figures"),
                                      trace_terms = c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a", "loglik"),
                                      make_traceplots_flag = TRUE) {
  make_dir(out_tables)
  make_dir(out_figures)
  
  message2("Building tables and figures from root: %s", root)
  
  missingness_tab <- make_missingness_table(root, out_tables)
  complete_case_tab <- make_complete_case_table(root, out_tables)
  param_tab <- make_param_recovery_table(root, out_tables)
  selection_tab <- make_model_selection_table(root, out_tables)
  evidence_tab <- make_evidence_table(root, out_tables)
  fp_tab <- make_fixed_point_table(root, out_tables)
  mcmc_tab <- make_mcmc_diagnostics_table(root, out_tables)
  runtime_tab <- make_runtime_table(root, out_tables)
  
  save_plot_outputs(plot_model_selection_rates(selection_tab), "fig_model_selection_rates", out_figures, width = 11, height = 7)
  save_plot_outputs(plot_parameter_bias(param_tab), "fig_parameter_bias_M2", out_figures, width = 10, height = 5)
  save_plot_outputs(plot_parameter_rmse(param_tab), "fig_parameter_rmse_M2", out_figures, width = 10, height = 5)
  save_plot_outputs(plot_evidence_heatmap(evidence_tab), "fig_evidence_heatmap_vs_M2", out_figures, width = 10, height = 5.5)
  save_plot_outputs(plot_fixed_point_sensitivity(fp_tab), "fig_fixed_point_sensitivity", out_figures, width = 9, height = 5)
  save_plot_outputs(plot_mcmc_ess(mcmc_tab), "fig_mcmc_ess", out_figures, width = 9, height = 5)
  save_plot_outputs(plot_runtime(runtime_tab), "fig_runtime_by_workflow", out_figures, width = 10, height = 5.5)
  
  if (isTRUE(make_traceplots_flag)) {
    make_traceplots(root, out_figures, trace_terms = trace_terms)
  }
  
  manifest <- tibble::tibble(
    output_type = c("tables", "figures"),
    path = c(normalizePath(out_tables, winslash = "/", mustWork = FALSE),
             normalizePath(out_figures, winslash = "/", mustWork = FALSE))
  )
  readr::write_csv(manifest, file.path(root, "paper_outputs_manifest.csv"))
  
  message2("Done. Tables: %s", out_tables)
  message2("Done. Figures: %s", out_figures)
  
  invisible(list(
    missingness = missingness_tab,
    complete_case = complete_case_tab,
    parameter_recovery = param_tab,
    model_selection = selection_tab,
    evidence = evidence_tab,
    fixed_point = fp_tab,
    mcmc = mcmc_tab,
    runtime = runtime_tab,
    out_tables = out_tables,
    out_figures = out_figures
  ))
}

# Uncomment for interactive use:
# make_paper_tables_figures(root = "outputs/smoke10")
