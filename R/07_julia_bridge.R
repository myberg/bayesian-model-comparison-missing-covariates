# -----------------------------------------------------------------------------
# 07_julia_bridge.R
# Bridge between the R simulation/ML workflow and the Julia Bayesian workflow.
# -----------------------------------------------------------------------------

# Export simulated datasets generated in R so that Julia can work on exactly the
# same data as the ML / MICE / missForest benchmarks.
export_simulation_data <- function(n_rep,
                                            mechanisms = c("MCAR", "MAR", "MNAR"),
                                            sim_control_arg = NULL,
                                            model_specs_arg = NULL,
                                            out_dir = file.path("data", "simulated")) {
  if (is.null(sim_control_arg)) {
    if (exists("sim_control", envir = .GlobalEnv)) {
      sim_control_arg <- get("sim_control", envir = .GlobalEnv)
    } else {
      sim_control_arg <- NULL
    }
  }
  
  if (is.null(model_specs_arg)) {
    if (exists("model_specs", envir = .GlobalEnv)) {
      model_specs_arg <- get("model_specs", envir = .GlobalEnv)
    } else {
      stop("No model_specs supplied and no global object 'model_specs' found.")
    }
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  sim_data <- make_simulation_data(
    n_rep = n_rep,
    mechanisms = mechanisms,
    sim_control = sim_control_arg
  )
  
  index <- data.frame()
  
  for (i in seq_along(sim_data)) {
    rep_id <- sim_data[[i]]$rep_id
    mech <- sim_data[[i]]$mechanism
    
    file_full <- file.path(out_dir, sprintf("rep_%03d_%s_full.csv", rep_id, mech))
    file_mis  <- file.path(out_dir, sprintf("rep_%03d_%s_mis.csv",  rep_id, mech))
    
    utils::write.csv(sim_data[[i]]$data_full, file_full, row.names = FALSE, na = "")
    utils::write.csv(sim_data[[i]]$data_mis,  file_mis,  row.names = FALSE, na = "")
    
    index <- rbind(index, data.frame(
      rep_id = rep_id,
      mechanism = mech,
      file_full = normalizePath(file_full, winslash = "/", mustWork = FALSE),
      file_mis = normalizePath(file_mis, winslash = "/", mustWork = FALSE),
      complete_case_rate = sim_data[[i]]$complete_case_rate,
      stringsAsFactors = FALSE
    ))
  }
  
  index_file <- file.path(out_dir, "simulation_index.csv")
  utils::write.csv(index, index_file, row.names = FALSE)
  
  model_tab <- do.call(rbind, lapply(names(model_specs_arg), function(nm) {
    fixed_formula <- model_specs_arg[[nm]]$fixed
    terms <- attr(stats::terms(fixed_formula), "term.labels")
    
    data.frame(
      model = nm,
      terms = paste(terms, collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
  
  model_file <- file.path(out_dir, "model_specs.csv")
  utils::write.csv(model_tab, model_file, row.names = FALSE)
  
  if (exists("sim_config", envir = .GlobalEnv) && !is.null(sim_config$beta_true)) {
    true_beta <- sim_config$beta_true
  } else {
    true_beta <- c(
      "(Intercept)" = 0.20,
      "x2" = -1.50,
      "x3" = 0.80,
      "x4" = -0.30,
      "x5B" = -0.40,
      "x5C" = -0.30,
      "x7" = 0.00,
      "x8B" = 0.00,
      "x8C" = 0.00,
      "x10" = 0.00
    )
  }
  
  true_beta_file <- file.path(out_dir, "true_beta.csv")
  utils::write.csv(
    data.frame(term = names(true_beta), true_value = as.numeric(true_beta)),
    true_beta_file,
    row.names = FALSE
  )
  
  message("Exported ", nrow(index), " simulation datasets to ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  
  invisible(list(
    sim_data = sim_data,
    index = index,
    index_file = normalizePath(index_file, winslash = "/", mustWork = FALSE),
    model_file = normalizePath(model_file, winslash = "/", mustWork = FALSE),
    true_beta_file = normalizePath(true_beta_file, winslash = "/", mustWork = FALSE)
  ))
}

inspect_simulation_index <- function(index_file = file.path("data", "simulated", "simulation_index.csv"),
                                     expected_n_rep = NULL,
                                     expected_mechanisms = NULL,
                                     stop_on_mismatch = FALSE) {
  if (!file.exists(index_file)) {
    stop("Simulation index not found: ", index_file,
         ". Run export_simulation_data() first.")
  }
  
  idx <- utils::read.csv(index_file, stringsAsFactors = FALSE)
  msg <- paste0(
    "Simulation index contains ", nrow(idx), " dataset rows; reps = ",
    paste(sort(unique(idx$rep_id)), collapse = ","),
    "; mechanisms = ", paste(sort(unique(idx$mechanism)), collapse = ",")
  )
  message(msg)
  
  if (!is.null(expected_n_rep) || !is.null(expected_mechanisms)) {
    exp_rep <- if (is.null(expected_n_rep)) length(unique(idx$rep_id)) else as.integer(expected_n_rep)
    exp_mech <- if (is.null(expected_mechanisms)) unique(idx$mechanism) else expected_mechanisms
    expected_rows <- exp_rep * length(exp_mech)
    ok <- nrow(idx) >= expected_rows &&
      all(seq_len(exp_rep) %in% unique(idx$rep_id)) &&
      all(exp_mech %in% unique(idx$mechanism))
    
    if (!ok) {
      warn <- paste0(
        "The simulation index does not match the requested design. Expected at least ",
        expected_rows, " rows (= ", exp_rep, " replications x ",
        length(exp_mech), " mechanisms), but found ", nrow(idx),
        ". Re-run export_simulation_data(n_rep = ..., mechanisms = ...)."
      )
      if (isTRUE(stop_on_mismatch)) stop(warn) else warning(warn)
    }
  }
  
  invisible(idx)
}

run_julia_bayes <- function(
  index_file = file.path("data", "simulated", "simulation_index.csv"),
  model_file = file.path("data", "simulated", "model_specs.csv"),
  julia_project = "julia",
  julia_script = file.path("julia", "run_bayes.jl"),
  out_dir = file.path("output", "julia_bayes"),
  julia_bin = Sys.getenv("JULIA_BIN", unset = "julia"),
  n_iter = 4000L,
  burnin = 1000L,
  thin = 5L,
  tau2_beta = 1.0,
  c0 = 2.0,
  d0 = 1.0,
  impute_mode = c("mean_once", "cart_once", "cart_every", "cartbb_once", "cartbb_every", "none"),
  data_version = c("mis", "full", "cc"),
  cart_every = 25L,
  cart_max_depth = 5L,
  cart_min_samples_leaf = 10L,
  cart_max_donors = 50L,
  conditional_medoid = FALSE,
  medoid_reference_model = "M2_true",
  medoid_n_iter = n_iter,
  medoid_thin = thin,
  chib_reduced_iter = 2000L,
  chib_reduced_burnin = 500L,
  chib_reduced_thin = 2L,
  compute_evidence = TRUE,
  verbose_every = 500L,
  models = NULL,
  reps = NULL,
  mechanisms_filter = NULL,
  expected_n_rep = NULL,
  expected_mechanisms = NULL,
  check_index = TRUE,
  export_draws = TRUE,
  diagnostics_dir = NULL,
  draw_reps = 1L,
  draw_models = "M2_true",
  draw_mechanisms = c("MCAR", "MAR"),
  extra_args = character()) {
  impute_mode <- match.arg(impute_mode)
  data_version <- match.arg(data_version)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  export_draws <- isTRUE(export_draws) || tolower(as.character(export_draws)) %in% c("true", "1", "yes")
  
  if (is.null(diagnostics_dir)) {
    out_norm <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
    
    if (grepl("/99_raw/", out_norm, fixed = TRUE)) {
      root_dir <- sub("/99_raw/.*$", "", out_norm)
    } else {
      root_dir <- dirname(out_norm)
    }
    
    diagnostics_dir <- file.path(root_dir, "06_mcmc_diagnostics")
  }
  
  diagnostics_dir <- normalizePath(
    diagnostics_dir,
    winslash = "/",
    mustWork = FALSE
  )
  
  
  if (isTRUE(check_index)) {
    inspect_simulation_index(
      index_file = index_file,
      expected_n_rep = expected_n_rep,
      expected_mechanisms = expected_mechanisms,
      stop_on_mismatch = !is.null(expected_n_rep) || !is.null(expected_mechanisms)
    )
  }
  
  quote_path_arg <- function(x) {
    x <- gsub("\\\\", "/", x)
    shQuote(x, type = if (.Platform$OS.type == "windows") "cmd" else "sh")
  }
  
  args <- c(
    paste0("--project=", julia_project),
    julia_script,
    
    "--index", quote_path_arg(index_file),
    "--models", quote_path_arg(model_file),
    "--outdir", quote_path_arg(out_dir),
    
    "--n_iter", as.character(as.integer(n_iter)),
    "--burnin", as.character(as.integer(burnin)),
    "--thin", as.character(as.integer(thin)),
    "--tau2_beta", as.character(tau2_beta),
    "--c0", as.character(c0),
    "--d0", as.character(d0),
    "--impute_mode", impute_mode,
    "--data_version", data_version,
    "--cart_every", as.character(as.integer(cart_every)),
    "--cart_max_depth", as.character(as.integer(cart_max_depth)),
    "--cart_min_samples_leaf", as.character(as.integer(cart_min_samples_leaf)),
    "--cart_max_donors", as.character(as.integer(cart_max_donors)),
    "--conditional_medoid", tolower(as.character(isTRUE(conditional_medoid))),
    "--medoid_reference_model", medoid_reference_model,
    "--medoid_n_iter", as.character(as.integer(medoid_n_iter)),
    "--medoid_burnin", as.character(as.integer(medoid_burnin)),
    "--medoid_thin", as.character(as.integer(medoid_thin)),
    "--chib_reduced_iter", as.character(as.integer(chib_reduced_iter)),
    "--chib_reduced_burnin", as.character(as.integer(chib_reduced_burnin)),
    "--chib_reduced_thin", as.character(as.integer(chib_reduced_thin)),
    "--compute_evidence", tolower(as.character(isTRUE(compute_evidence))),
    "--verbose_every", as.character(as.integer(verbose_every)),
    
    "--export_draws", tolower(as.character(export_draws)),
    "--diagnostics_dir", quote_path_arg(diagnostics_dir),
    "--draw_reps", paste(draw_reps, collapse = ","),
    "--draw_models", paste(draw_models, collapse = ","),
    "--draw_mechanisms", paste(draw_mechanisms, collapse = ",")
  )
  
  if (!is.null(models) && length(models) > 0L) {
    args <- c(args, "--only_models", paste(models, collapse = ","))
  }
  
  if (!is.null(reps) && length(reps) > 0L) {
    args <- c(args, "--only_reps", paste(as.integer(reps), collapse = ","))
  }
  
  if (!is.null(mechanisms_filter) && length(mechanisms_filter) > 0L) {
    args <- c(args, "--only_mechanisms", paste(mechanisms_filter, collapse = ","))
  }
  
  args <- c(args, extra_args)
  
  message("Running Julia: ", julia_bin, " ", paste(args, collapse = " "))
  status <- system2(command = julia_bin, args = args)
  
  if (!identical(status, 0L)) {
    stop("Julia Bayesian workflow failed with exit status ", status, ".")
  }
  
  invisible(file.path(out_dir, "bayes_model_results.csv"))
}

read_julia_bayes_results <- function(out_dir = file.path("output", "julia_bayes")) {
  model_file <- file.path(out_dir, "bayes_model_results.csv")
  estimate_file <- file.path(out_dir, "bayes_estimates.csv")
  marglik_file <- file.path(out_dir, "bayes_marginal_likelihood.csv")
  predictive_file <- file.path(out_dir, "bayes_predictive_criteria.csv")
  medoid_metadata_file <- file.path(out_dir, "medoid_imputation_metadata.csv")
  
  safe_read_csv <- function(path) {
    if (!file.exists(path)) return(NULL)
    info <- file.info(path)
    if (is.na(info$size) || info$size == 0L) return(NULL)
    
    txt <- readLines(path, warn = FALSE)
    if (length(txt) == 0L) return(NULL)
    if (length(txt) == 1L && !nzchar(trimws(txt[1L]))) return(NULL)
    
    tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE),
      error = function(e) {
        warning("Could not read Julia CSV file '", path, "': ", conditionMessage(e))
        NULL
      }
    )
  }
  
  list(
    bayes = safe_read_csv(model_file),
    bayes_estimates = safe_read_csv(estimate_file),
    bayes_marginal_likelihood = safe_read_csv(marglik_file),
    bayes_predictive_criteria = safe_read_csv(predictive_file),
    medoid_imputation_metadata = safe_read_csv(medoid_metadata_file)
  )
}

summarise_predictive_criteria <- function(criteria) {
  if (is.null(criteria) || nrow(criteria) == 0L) return(NULL)
  
  num_vars <- c("lppd", "p_waic", "elpd_waic", "WAIC", "elpd_loo", "p_loo", "LOO")
  num_vars <- num_vars[num_vars %in% names(criteria)]
  group_vars <- c("method", "model", "mechanism", "data_version", "impute_mode")
  group_vars <- group_vars[group_vars %in% names(criteria)]
  
  split_fac <- interaction(criteria[group_vars], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(criteria, split_fac), function(x) {
    g <- x[1L, group_vars, drop = FALSE]
    vals <- list(g)
    for (v in num_vars) {
      vals[[paste0("mean_", v)]] <- mean(x[[v]], na.rm = TRUE)
      vals[[paste0("sd_", v)]] <- stats::sd(x[[v]], na.rm = TRUE)
      vals[[paste0("median_", v)]] <- stats::median(x[[v]], na.rm = TRUE)
    }
    vals[["n_rep"]] <- length(unique(x$rep))
    as.data.frame(vals, stringsAsFactors = FALSE)
  })
  
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

make_bayes_selection_table <- function(model_results,
                                       criterion = c("log_marginal_likelihood", "WAIC", "LOO", "logLik_mean"),
                                       reference_model = "M2_true") {
  criterion <- match.arg(criterion)
  if (is.null(model_results) || nrow(model_results) == 0L) return(NULL)
  if (!criterion %in% names(model_results)) {
    stop("Criterion '", criterion, "' not found in model_results.")
  }
  
  x <- model_results[!is.na(model_results[[criterion]]), , drop = FALSE]
  if (nrow(x) == 0L) return(NULL)
  
  lower_is_better <- criterion %in% c("WAIC", "LOO", "AIC", "BIC")
  group_vars <- c("method", "rep", "mechanism", "data_version", "impute_mode")
  group_vars <- group_vars[group_vars %in% names(x)]
  split_fac <- interaction(x[group_vars], drop = TRUE, lex.order = TRUE)
  
  out <- lapply(split(x, split_fac), function(z) {
    if (nrow(z) == 0L) return(NULL)
    value <- z[[criterion]]
    idx <- if (lower_is_better) which.min(value) else which.max(value)
    selected <- z[idx, , drop = FALSE]
    selected$criterion <- criterion
    selected$selected_model <- selected$model
    selected$selected_value <- selected[[criterion]]
    selected$reference_model <- reference_model
    selected$selected_reference <- selected$model == reference_model
    selected
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

summarise_bayes_selection <- function(selection_tab) {
  if (is.null(selection_tab) || nrow(selection_tab) == 0L) return(NULL)
  
  group_vars <- c("method", "criterion", "mechanism", "data_version", "impute_mode", "selected_model")
  group_vars <- group_vars[group_vars %in% names(selection_tab)]
  
  counts <- as.data.frame(table(selection_tab[group_vars]), stringsAsFactors = FALSE)
  names(counts)[names(counts) == "Freq"] <- "n_selected"
  counts <- counts[counts$n_selected > 0L, , drop = FALSE]
  
  total_vars <- setdiff(group_vars, "selected_model")
  totals <- stats::aggregate(n_selected ~ ., data = counts[, c(total_vars, "n_selected"), drop = FALSE], sum)
  names(totals)[names(totals) == "n_selected"] <- "n_total"
  
  out <- merge(counts, totals, by = total_vars, all.x = TRUE)
  out$selection_rate <- out$n_selected / out$n_total
  out[order(out$mechanism, out$criterion, -out$selection_rate), , drop = FALSE]
}

summarise_bayes_marginal_likelihood <- function(marglik, fixed_point = "mean") {
  if (is.null(marglik) || nrow(marglik) == 0L) return(NULL)
  x <- marglik[marglik$fixed_point == fixed_point, , drop = FALSE]
  if (nrow(x) == 0L) return(NULL)
  
  stats::aggregate(
    log_marginal_likelihood ~ method + model + mechanism + data_version + impute_mode + fixed_point,
    data = x,
    FUN = function(z) c(mean = mean(z, na.rm = TRUE),
                        sd = stats::sd(z, na.rm = TRUE),
                        median = stats::median(z, na.rm = TRUE),
                        n_rep = sum(!is.na(z)))
  ) |>
    as.data.frame() |>
    within({
      mean_log_marginal_likelihood <- log_marginal_likelihood[, "mean"]
      sd_log_marginal_likelihood <- log_marginal_likelihood[, "sd"]
      median_log_marginal_likelihood <- log_marginal_likelihood[, "median"]
      n_rep <- log_marginal_likelihood[, "n_rep"]
      log_marginal_likelihood <- NULL
    })
}

make_bayes_factor_table <- function(marglik, reference_model = "M2_true", fixed_point = "mean") {
  if (is.null(marglik) || nrow(marglik) == 0L) return(NULL)
  x <- marglik[marglik$fixed_point == fixed_point, , drop = FALSE]
  if (nrow(x) == 0L) return(NULL)
  
  ref <- x[x$model == reference_model,
           c("method", "rep", "mechanism", "data_version", "impute_mode",
             "log_marginal_likelihood"),
           drop = FALSE]
  names(ref)[names(ref) == "log_marginal_likelihood"] <- "log_marginal_likelihood_ref"
  
  out <- merge(
    x, ref,
    by = c("method", "rep", "mechanism", "data_version", "impute_mode"),
    all.x = TRUE
  )
  out$reference_model <- reference_model
  out$lnBF_model_vs_ref <- out$log_marginal_likelihood - out$log_marginal_likelihood_ref
  out$delta_log_marginal_likelihood_ref_minus_model <-
    out$log_marginal_likelihood_ref - out$log_marginal_likelihood
  out
}

summarise_bayes_factor_table <- function(bf_tab) {
  if (is.null(bf_tab) || nrow(bf_tab) == 0L) return(NULL)
  stats::aggregate(
    lnBF_model_vs_ref ~ method + model + reference_model + mechanism + data_version + impute_mode + fixed_point,
    data = bf_tab,
    FUN = function(z) c(mean = mean(z, na.rm = TRUE),
                        sd = stats::sd(z, na.rm = TRUE),
                        median = stats::median(z, na.rm = TRUE),
                        n_rep = sum(!is.na(z)))
  ) |>
    as.data.frame() |>
    within({
      mean_lnBF_model_vs_ref <- lnBF_model_vs_ref[, "mean"]
      sd_lnBF_model_vs_ref <- lnBF_model_vs_ref[, "sd"]
      median_lnBF_model_vs_ref <- lnBF_model_vs_ref[, "median"]
      n_rep <- lnBF_model_vs_ref[, "n_rep"]
      lnBF_model_vs_ref <- NULL
    })
}

combine_model_results <- function(ml_results, julia_results) {
  add_missing_cols <- function(x, cols) {
    missing_cols <- setdiff(cols, names(x))
    for (nm in missing_cols) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  }
  
  if (is.null(ml_results)) return(julia_results)
  if (is.null(julia_results)) return(ml_results)
  
  all_cols <- union(names(ml_results), names(julia_results))
  rbind(add_missing_cols(ml_results, all_cols), add_missing_cols(julia_results, all_cols))
}


# Harmonise ML and Julia/Bayes estimates and summarise across Monte Carlo replications.
harmonise_estimates_for_summary <- function(ml_estimates = NULL,
                                            bayes_estimates = NULL,
                                            true_values = NULL) {
  if (is.null(true_values)) {
    true_values <- c(
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
  
  out <- list()
  
  if (!is.null(ml_estimates) && nrow(ml_estimates) > 0L) {
    ml <- ml_estimates
    ml$estimate_common <- ml$estimate
    ml$se_common <- ml$std_error
    ml$source <- "ML"
    if (!"impute_mode" %in% names(ml)) ml$impute_mode <- NA_character_
    if (!"data_version" %in% names(ml)) ml$data_version <- NA_character_
    out[[length(out) + 1L]] <- ml
  }
  
  if (!is.null(bayes_estimates) && nrow(bayes_estimates) > 0L) {
    by <- bayes_estimates
    by$estimate_common <- by$posterior_mean
    by$se_common <- by$posterior_sd
    by$source <- "Bayes"
    if (!"error" %in% names(by)) by$error <- NA_character_
    if (!"std_error" %in% names(by)) by$std_error <- NA_real_
    if (!"z_value" %in% names(by)) by$z_value <- NA_real_
    if (!"p_value" %in% names(by)) by$p_value <- NA_real_
    out[[length(out) + 1L]] <- by
  }
  
  if (length(out) == 0L) return(NULL)
  
  all_cols <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(x) {
    miss <- setdiff(all_cols, names(x))
    for (cc in miss) x[[cc]] <- NA
    x[, all_cols, drop = FALSE]
  })
  
  dat <- do.call(rbind, out)
  dat$true_value <- unname(true_values[dat$term])
  dat$bias <- dat$estimate_common - dat$true_value
  dat$abs_bias <- abs(dat$bias)
  dat$sq_error <- dat$bias^2
  dat
}

summarise_estimates_over_reps <- function(estimates_long,
                                          group_vars = c("source", "method", "model", "mechanism",
                                                         "data_version", "impute_mode", "term", "component")) {
  if (is.null(estimates_long) || nrow(estimates_long) == 0L) return(NULL)
  
  group_vars <- group_vars[group_vars %in% names(estimates_long)]
  split_fac <- interaction(estimates_long[group_vars], drop = TRUE, lex.order = TRUE)
  
  pieces <- lapply(split(estimates_long, split_fac), function(x) {
    g <- x[1L, group_vars, drop = FALSE]
    data.frame(
      g,
      true_value = { tv <- x$true_value[!is.na(x$true_value)]; if (length(tv) == 0L) NA_real_ else tv[1L] },
      mean_estimate = mean(x$estimate_common, na.rm = TRUE),
      mc_sd = stats::sd(x$estimate_common, na.rm = TRUE),
      mean_se = mean(x$se_common, na.rm = TRUE),
      bias = mean(x$bias, na.rm = TRUE),
      mean_abs_bias = mean(x$abs_bias, na.rm = TRUE),
      rmse = sqrt(mean(x$sq_error, na.rm = TRUE)),
      n_rep = length(unique(x$rep)),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}
