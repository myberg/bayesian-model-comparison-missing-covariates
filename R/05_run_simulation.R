# -----------------------------------------------------------------------------
# 05_run_simulation.R
# Driver script for the additional reviewer-requested analyses
# -----------------------------------------------------------------------------

run_ml_benchmarks_one <- function(
  sim_item,
  model_specs,
  run_missforest = TRUE,
  seed = 123L,
  id_var = "id",
  y_var = "y",
  mice_m = 5L,
  mice_maxit = 5L) {
  
  if (is.null(model_specs) || length(model_specs) == 0L) {
    stop("run_ml_benchmarks_one: model_specs is empty or NULL.")
  }
  
  if (is.null(names(model_specs)) || any(names(model_specs) == "")) {
    stop("run_ml_benchmarks_one: model_specs must be a named list.")
  }
  
  if (!"data_mis" %in% names(sim_item)) {
    stop(
      "run_ml_benchmarks_one: sim_item has no element 'data_mis'. Available elements are: ",
      paste(names(sim_item), collapse = ", ")
    )
  }
  
  rep_id <- sim_item$rep_id
  mech <- sim_item$mechanism
  
  data_mis <- as.data.frame(sim_item$data_mis)
  
  if (!y_var %in% names(data_mis)) {
    stop(
      "run_ml_benchmarks_one: y_var '", y_var,
      "' not found in data_mis. Available columns: ",
      paste(names(data_mis), collapse = ", ")
    )
  }
  
  if (!id_var %in% names(data_mis)) {
    stop(
      "run_ml_benchmarks_one: id_var '", id_var,
      "' not found in data_mis. Available columns: ",
      paste(names(data_mis), collapse = ", ")
    )
  }
  
  data_mis[[id_var]] <- as.factor(data_mis[[id_var]])
  
  all_model_vars <- unique(unlist(lapply(model_specs, function(x) {
    all.vars(x$fixed)
  })))
  
  common_cc_vars <- unique(c(y_var, id_var, all_model_vars))
  common_cc_vars <- common_cc_vars[common_cc_vars %in% names(data_mis)]
  
  cc_rows_common <- stats::complete.cases(data_mis[, common_cc_vars, drop = FALSE])
  data_cc_common <- data_mis[cc_rows_common, , drop = FALSE]
  
  message(
    "Common complete-case sample size = ",
    nrow(data_cc_common), " / ", nrow(data_mis)
  )
  
  # setup counters
  estimate_results <- list()
  estimate_counter <- 1L
  results <- list()
  counter <- 1L
  
  
  # for loop over all models
  for (model_name in names(model_specs)) {
    
    fixed_formula <- model_specs[[model_name]]$fixed
    
    message("ML benchmark | model = ", model_name)
    
    # Complete cases
    cc <- fit_complete_case_ml(
      data = data_cc_common,
      fixed_formula = fixed_formula,
      id_var = id_var,
      y_var = y_var
    )
    
    results[[counter]] <- extract_ml_fit(
      fit = cc$value,
      method = "Complete case + ML",
      model_name = model_name,
      rep = rep_id,
      mechanism = mech,
      runtime_seconds = cc$runtime_seconds
    )
    counter <- counter + 1L
    
    estimate_results[[estimate_counter]] <- extract_glmer_estimates(
      fit = cc$value,
      method = "Complete case + ML",
      model = model_name,
      rep = rep_id,
      mechanism = mech,
      id_var = id_var
    )
    estimate_counter <- estimate_counter + 1L
    
    
    # MICE + ML
    mice_res <- tryCatch(
      fit_mice_ml(
        data = data_mis,
        fixed_formula = fixed_formula,
        id_var = id_var,
        y_var = y_var,
        m = mice_m,
        maxit = mice_maxit,
        seed = seed + rep_id
      ),
      error = function(e) e
    )
    
    if (inherits(mice_res, "error")) {
      results[[counter]] <- data.frame(
        method = "MICE-CART + ML",
        model = model_name,
        rep = rep_id,
        mechanism = mech,
        converged = FALSE,
        error = mice_res$message,
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        nobs = NA_integer_,
        runtime_seconds = NA_real_
      )
    } else {
      first_fit <- mice_res$fits[[1]]
      
      results[[counter]] <- data.frame(
        method = "MICE-CART + ML",
        model = model_name,
        rep = rep_id,
        mechanism = mech,
        converged = TRUE,
        error = NA_character_,
        logLik = mice_res$logLik,
        AIC = mice_res$AIC,
        BIC = mice_res$BIC,
        nobs = if (!is.null(first_fit)) stats::nobs(first_fit) else NA_integer_,
        runtime_seconds = mice_res$runtime_seconds
      )
    }
    
    counter <- counter + 1L
    
    estimate_results[[estimate_counter]] <- tryCatch(
      extract_glmer_estimates(
        fit = cc$value,
        method = "Complete case + ML",
        model = model_name,
        rep = rep_id,
        mechanism = mech,
        id_var = id_var
      ),
      error = function(e) {
        data.frame(
          method = "Complete case + ML",
          model = model_name,
          rep = rep_id,
          mechanism = mech,
          term = NA_character_,
          estimate = NA_real_,
          std_error = NA_real_,
          z_value = NA_real_,
          p_value = NA_real_,
          component = "fixed",
          error = paste("Estimate extraction failed:", e$message),
          stringsAsFactors = FALSE
        )
      }
    )
    
    estimate_counter <- estimate_counter + 1L
    
    
    # missForest + ML
    if (isTRUE(run_missforest)) {
      
      mf <- tryCatch(
        fit_missforest_ml(
          data = data_mis,
          fixed_formula = fixed_formula,
          id_var = id_var,
          y_var = y_var,
          seed = seed + 1000L + rep_id
        ),
        error = function(e) list(value = e, runtime_seconds = NA_real_)
      )
      
      mf_fit <- if (inherits(mf$value, "error")) {
        mf$value
      } else {
        mf$value$fit
      }
      
      results[[counter]] <- extract_ml_fit(
        fit = mf_fit,
        method = "missForest + ML",
        model_name = model_name,
        rep = rep_id,
        mechanism = mech,
        runtime_seconds = mf$runtime_seconds
      )
      
      counter <- counter + 1L
    }
  }
  
  if (length(results) == 0L) {
    stop("No ML results were produced. Check whether model_specs is empty.")
  }
  
  list(
    metrics = do.call(rbind, results),
    estimates = if (length(estimate_results) > 0L) {
      do.call(rbind, estimate_results)
    } else {
      NULL
    }
  )
}


run_bayes_one <- function(
  sim_item,
  model_specs,
  id_var = "id",
  y_var = "y",
  mcmc_control,
  tau2_beta = 1,
  seed = 123L
) {
  
  if (is.null(model_specs) || length(model_specs) == 0L) {
    stop("run_bayes_one: model_specs is empty or NULL.")
  }
  
  data_mis <- as.data.frame(sim_item$data_mis)
  rep_id <- sim_item$rep_id
  mech <- sim_item$mechanism
  
  results <- list()
  fits <- list()
  
  for (model_name in names(model_specs)) {
    
    fixed_formula <- model_specs[[model_name]]$fixed
    tmp_design <- model_matrix_from_formula(fixed_formula, data_mis)
    prior <- make_beta_prior(ncol(tmp_design), tau2_beta = tau2_beta)
    
    timed <- with_runtime({
      fit_bayes_re_probit_cart(
        data = data_mis,
        fixed_formula = fixed_formula,
        id_var = id_var,
        y_var = y_var,
        n_iter = mcmc_control$n_iter,
        burnin = mcmc_control$burnin,
        beta_prior = prior,
        c0 = sigma2_a_prior$c0,
        d0 = sigma2_a_prior$d0,
        cart_control = list(
          minbucket_range = mcmc_control$cart_minbucket_range,
          maxdepth = mcmc_control$cart_maxdepth,
          cp = mcmc_control$cart_cp,
          max_donors = mcmc_control$cart_max_donors
        ),
        seed = seed + 2000L + rep_id,
        verbose_every = mcmc_control$verbose_every
      )
    })
    
    fits[[model_name]] <- timed$value
    diag <- summarise_mcmc_diagnostics(timed$value)
    post <- posterior_draws_after_burnin(timed$value)
    
    results[[model_name]] <- data.frame(
      method = "Bayesian CART-assisted MCMC",
      model = model_name,
      rep = rep_id,
      mechanism = mech,
      tau2_beta = tau2_beta,
      logLik_mean = mean(post$log_lik),
      logLik_median = median(post$log_lik),
      min_ESS = min(diag$ess, na.rm = TRUE),
      runtime_seconds = timed$runtime_seconds
    )
  }
  
  list(
    summary = do.call(rbind, results),
    fits = fits
  )
}


run_simulation <- function(
  n_rep = 1L,
  mechanisms = c("MCAR", "MAR", "MNAR"),
  run_ml = TRUE,
  run_bayes = FALSE,
  run_missforest = FALSE,
  sim_control_arg = NULL,
  model_specs_arg = NULL,
  seed = 202605L
) {
  
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
  
  if (is.null(model_specs_arg) || length(model_specs_arg) == 0L) {
    stop("model_specs_arg is empty or NULL.")
  }
  
  sim_data <- make_simulation_data(
    n_rep = n_rep,
    mechanisms = mechanisms,
    sim_control = sim_control_arg
  )
  
  ml_results <- list()
  ml_estimates <- list()
  bayes_results <- list()
  bayes_fits <- list()
  
  for (i in seq_along(sim_data)) {
    
    message(
      "Simulation item ", i, " / ", length(sim_data),
      " | rep=", sim_data[[i]]$rep_id,
      " | mechanism=", sim_data[[i]]$mechanism
    )
    
    if (isTRUE(run_ml)) {
      ml_tmp <- run_ml_benchmarks_one(
        sim_item = sim_data[[i]],
        model_specs = model_specs_arg,
        run_missforest = run_missforest,
        seed = seed
      )
      
      ml_results[[i]] <- ml_tmp$metrics
      ml_estimates[[i]] <- ml_tmp$estimates
      
      
    }
    
    if (isTRUE(run_bayes)) {
      br <- run_bayes_one(
        sim_item = sim_data[[i]],
        model_specs = model_specs_arg,
        mcmc_control = mcmc_control,
        tau2_beta = 1,
        seed = seed
      )
      
      bayes_results[[i]] <- br$summary
      bayes_fits[[i]] <- br$fits
    }
  }
  
  # output
  out <- list(
    sim_data = sim_data,
    ml = if (length(ml_results) > 0L) do.call(rbind, ml_results) else NULL,
    ml_estimates = if (length(ml_estimates) > 0L) do.call(rbind, ml_estimates) else NULL,
    bayes = if (length(bayes_results) > 0L) do.call(rbind, bayes_results) else NULL,
    bayes_fits = bayes_fits
  )
  
  if (exists("output_dir", envir = .GlobalEnv)) {
    out_dir <- file.path(output_dir, "rds")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    saveRDS(out, file.path(out_dir, "simulation_results.rds"))
  }
  
  out
}