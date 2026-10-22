# -----------------------------------------------------------------------------
# 06_collect_tables.R
# Summaries for manuscript tables
# -----------------------------------------------------------------------------


summarise_ml_table <- function(ml_results, true_model = "M2_true") {
  selected <- select_model_by_min_bic(ml_results)
  if (is.null(selected) || nrow(selected) == 0) return(NULL)

  correct <- aggregate(model ~ method + mechanism, selected, function(x) mean(x == true_model, na.rm = TRUE))
  names(correct)[names(correct) == "model"] <- "correct_selection_rate"

  runt <- aggregate(runtime_seconds ~ method + mechanism, ml_results, function(x) mean(x, na.rm = TRUE))
  names(runt)[names(runt) == "runtime_seconds"] <- "mean_runtime_seconds"

  merge(correct, runt, by = c("method", "mechanism"), all = TRUE)
}

summarise_bayes_table <- function(bayes_results, true_model = "M2_true") {
  if (is.null(bayes_results) || nrow(bayes_results) == 0) return(NULL)

  aggregate(cbind(logLik_mean, min_ESS, runtime_seconds) ~ method + mechanism + model,
            bayes_results, function(x) mean(x, na.rm = TRUE))
}

make_bic_logbf_summary <- function(ml_results, reference_model = "M2_true") {
  x <- bic_logbf_table(ml_results, reference_model = reference_model)
  aggregate(cbind(delta_BIC_vs_ref, logBF_vs_ref) ~ method + mechanism + model,
            x, function(z) mean(z, na.rm = TRUE))
}

write_tables <- function(results,
                                  output_dir_arg = NULL) {
  
  if (is.null(output_dir_arg)) {
    if (exists("output_dir", envir = .GlobalEnv)) {
      output_dir_arg <- get("output_dir", envir = .GlobalEnv)
    } else {
      output_dir_arg <- "output"
    }
  }
  
  tables_dir <- file.path(output_dir_arg, "tables")
  
  if (!dir.exists(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE)
  }
  
  if (is.null(results$ml)) {
    warning("No ML results found in results$ml. No ML tables written.")
    return(invisible(NULL))
  }
  
  ml <- results$ml
  
  # Full ML results
  utils::write.csv(
    ml,
    file.path(tables_dir, "table_ml_full_results.csv"),
    row.names = FALSE
  )
  
  # Model selected by minimum BIC within method x mechanism x replication
  ml_ok <- ml[!is.na(ml$BIC), , drop = FALSE]
  
  if (nrow(ml_ok) > 0L) {
    split_ml <- split(
      ml_ok,
      list(ml_ok$method, ml_ok$mechanism, ml_ok$rep),
      drop = TRUE
    )
    
    selection_tab <- do.call(rbind, lapply(split_ml, function(d) {
      d[which.min(d$BIC), , drop = FALSE]
    }))
    
    selection_tab <- selection_tab[, c(
      "method", "mechanism", "rep", "model",
      "BIC", "AIC", "logLik", "nobs", "runtime_seconds",
      "converged", "error"
    ), drop = FALSE]
    
    utils::write.csv(
      selection_tab,
      file.path(tables_dir, "table_ml_selection_runtime.csv"),
      row.names = FALSE
    )
  } else {
    warning("No non-missing BIC values found. Selection table not written.")
    selection_tab <- NULL
  }
  
  # Average runtime by method and mechanism
  runtime_tab <- aggregate(
    runtime_seconds ~ method + mechanism,
    data = ml,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  names(runtime_tab)[names(runtime_tab) == "runtime_seconds"] <- "mean_runtime_seconds"
  
  utils::write.csv(
    runtime_tab,
    file.path(tables_dir, "table_runtime_summary.csv"),
    row.names = FALSE
  )
  
  invisible(list(
    ml_full = ml,
    selection = selection_tab,
    runtime = runtime_tab
  ))
}




extract_glmer_estimates <- function(fit,
                                    method,
                                    model,
                                    rep,
                                    mechanism,
                                    id_var = "id") {
  
  # Failed fit
  if (inherits(fit, "error")) {
    return(data.frame(
      method = method,
      model = model,
      rep = rep,
      mechanism = mechanism,
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      component = "fixed",
      error = fit$message,
      stringsAsFactors = FALSE
    ))
  }
  
  if (is.null(fit)) {
    return(data.frame(
      method = method,
      model = model,
      rep = rep,
      mechanism = mechanism,
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      component = "fixed",
      error = "Fit object is NULL.",
      stringsAsFactors = FALSE
    ))
  }
  
  # Prefer coefficient table from summary()
  coef_tab <- tryCatch(
    summary(fit)$coefficients,
    error = function(e) NULL
  )
  
  if (!is.null(coef_tab) && nrow(coef_tab) > 0L) {
    
    term <- rownames(coef_tab)
    
    estimate <- as.numeric(coef_tab[, "Estimate"])
    
    std_error <- if ("Std. Error" %in% colnames(coef_tab)) {
      as.numeric(coef_tab[, "Std. Error"])
    } else {
      rep(NA_real_, length(estimate))
    }
    
    z_value <- if ("z value" %in% colnames(coef_tab)) {
      as.numeric(coef_tab[, "z value"])
    } else {
      estimate / std_error
    }
    
    p_value <- if ("Pr(>|z|)" %in% colnames(coef_tab)) {
      as.numeric(coef_tab[, "Pr(>|z|)"])
    } else {
      2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
    }
    
  } else {
    
    beta_hat <- tryCatch(
      lme4::fixef(fit),
      error = function(e) numeric(0)
    )
    
    if (length(beta_hat) == 0L) {
      return(data.frame(
        method = method,
        model = model,
        rep = rep,
        mechanism = mechanism,
        term = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        z_value = NA_real_,
        p_value = NA_real_,
        component = "fixed",
        error = "No fixed-effect estimates could be extracted.",
        stringsAsFactors = FALSE
      ))
    }
    
    term <- names(beta_hat)
    if (is.null(term) || length(term) == 0L) {
      term <- paste0("beta_", seq_along(beta_hat))
    }
    
    estimate <- as.numeric(beta_hat)
    
    std_error <- tryCatch(
      {
        V <- as.matrix(stats::vcov(fit))
        se <- sqrt(diag(V))
        as.numeric(se[term])
      },
      error = function(e) rep(NA_real_, length(estimate))
    )
    
    if (length(std_error) != length(estimate)) {
      std_error <- rep(NA_real_, length(estimate))
    }
    
    z_value <- estimate / std_error
    p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  }
  
  # Fixed effects
  out_fixed <- data.frame(
    method = rep(method, length(term)),
    model = rep(model, length(term)),
    rep = rep(rep, length(term)),
    mechanism = rep(mechanism, length(term)),
    term = term,
    estimate = estimate,
    std_error = std_error,
    z_value = z_value,
    p_value = p_value,
    component = rep("fixed", length(term)),
    error = rep(NA_character_, length(term)),
    stringsAsFactors = FALSE
  )
  
  # Random-effect variance
  sigma2_a <- tryCatch(
    {
      vc <- lme4::VarCorr(fit)
      
      if (id_var %in% names(vc)) {
        as.numeric(vc[[id_var]][1, 1])
      } else if (length(vc) > 0L) {
        as.numeric(vc[[1]][1, 1])
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
  
  out_re <- data.frame(
    method = method,
    model = model,
    rep = rep,
    mechanism = mechanism,
    term = "sigma2_a",
    estimate = sigma2_a,
    std_error = NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    component = "random_effect_variance",
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  
  rbind(out_fixed, out_re)
}



pool_mice_glmer_estimates <- function(fits,
                                      method,
                                      model,
                                      rep,
                                      mechanism,
                                      id_var = "id") {
  fits_ok <- fits[!vapply(fits, inherits, logical(1), what = "error")]
  
  if (length(fits_ok) == 0L) {
    return(data.frame(
      method = method,
      model = model,
      rep = rep,
      mechanism = mechanism,
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      component = "fixed",
      error = "All MICE fits failed."
    ))
  }
  
  terms <- names(lme4::fixef(fits_ok[[1]]))
  m <- length(fits_ok)
  
  Q <- do.call(rbind, lapply(fits_ok, function(f) {
    lme4::fixef(f)[terms]
  }))
  
  U_list <- lapply(fits_ok, function(f) {
    as.matrix(stats::vcov(f))[terms, terms, drop = FALSE]
  })
  
  qbar <- colMeans(Q)
  
  ubar <- Reduce("+", U_list) / m
  
  if (m > 1L) {
    b <- stats::var(Q)
  } else {
    b <- matrix(0, nrow = length(terms), ncol = length(terms))
    rownames(b) <- colnames(b) <- terms
  }
  
  total_var <- ubar + (1 + 1 / m) * b
  se <- sqrt(diag(total_var))
  
  z <- qbar / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  
  out_fixed <- data.frame(
    method = method,
    model = model,
    rep = rep,
    mechanism = mechanism,
    term = terms,
    estimate = as.numeric(qbar),
    std_error = as.numeric(se),
    z_value = as.numeric(z),
    p_value = as.numeric(p),
    component = "fixed",
    error = NA_character_
  )
  
  sigma2_a <- vapply(fits_ok, function(f) {
    as.numeric(lme4::VarCorr(f)[[id_var]][1, 1])
  }, numeric(1))
  
  out_re <- data.frame(
    method = method,
    model = model,
    rep = rep,
    mechanism = mechanism,
    term = "sigma2_a",
    estimate = mean(sigma2_a, na.rm = TRUE),
    std_error = stats::sd(sigma2_a, na.rm = TRUE),
    z_value = NA_real_,
    p_value = NA_real_,
    component = "random_effect_variance",
    error = NA_character_
  )
  
  rbind(out_fixed, out_re)
}


make_ml_estimate_table <- function(results,
                                   model_specs,
                                   id_var = "id",
                                   y_var = "y",
                                   mice_m = 5L,
                                   mice_maxit = 5L,
                                   seed = 123L,
                                   output_dir_arg = NULL) {
  est_list <- list()
  counter <- 1L
  
  for (i in seq_along(results$sim_data)) {
    sim_item <- results$sim_data[[i]]
    
    rep_id <- sim_item$rep_id
    mech <- sim_item$mechanism
    data_mis <- as.data.frame(sim_item$data_mis)
    
    all_model_vars <- unique(unlist(lapply(model_specs, function(x) {
      all.vars(x$fixed)
    })))
    
    common_cc_vars <- unique(c(y_var, id_var, all_model_vars))
    common_cc_vars <- common_cc_vars[common_cc_vars %in% names(data_mis)]
    
    cc_rows_common <- stats::complete.cases(data_mis[, common_cc_vars, drop = FALSE])
    data_cc_common <- data_mis[cc_rows_common, , drop = FALSE]
    
    data_cc_common[[id_var]] <- as.factor(data_cc_common[[id_var]])
    data_mis[[id_var]] <- as.factor(data_mis[[id_var]])
    
    for (model_name in names(model_specs)) {
      fixed_formula <- model_specs[[model_name]]$fixed
      
      message("Estimate table | rep=", rep_id,
              " | mechanism=", mech,
              " | model=", model_name)
      
      # Complete case
      fit_cc <- fit_glmer_safe(
        data = data_cc_common,
        fixed_formula = fixed_formula,
        id_var = id_var,
        y_var = y_var
      )
      
      est_list[[counter]] <- extract_glmer_estimates(
        fit = fit_cc,
        method = "Complete case + ML",
        model = model_name,
        rep = rep_id,
        mechanism = mech,
        id_var = id_var
      )
      counter <- counter + 1L
      
      # MICE-CART
      fit_mice <- tryCatch(
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
      
      if (inherits(fit_mice, "error")) {
        est_list[[counter]] <- data.frame(
          method = "MICE-CART + ML",
          model = model_name,
          rep = rep_id,
          mechanism = mech,
          term = NA_character_,
          estimate = NA_real_,
          std_error = NA_real_,
          z_value = NA_real_,
          p_value = NA_real_,
          component = "fixed",
          error = fit_mice$message
        )
      } else {
        est_list[[counter]] <- pool_mice_glmer_estimates(
          fits = fit_mice$fits,
          method = "MICE-CART + ML",
          model = model_name,
          rep = rep_id,
          mechanism = mech,
          id_var = id_var
        )
      }
      
      counter <- counter + 1L
    }
  }
  
  est_tab <- do.call(rbind, est_list)
  
  if (is.null(output_dir_arg)) {
    if (exists("output_dir", envir = .GlobalEnv)) {
      output_dir_arg <- get("output_dir", envir = .GlobalEnv)
    } else {
      output_dir_arg <- "output"
    }
  }
  
  tables_dir <- file.path(output_dir_arg, "tables")
  if (!dir.exists(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE)
  }
  
  utils::write.csv(
    est_tab,
    file.path(tables_dir, "table_ml_estimates.csv"),
    row.names = FALSE
  )
  
  est_tab
}



summarise_estimates_over_reps <- function(est_tab, true_beta = NULL) {
  
  est_ok <- est_tab[
    !is.na(est_tab$estimate) &
      est_tab$component == "fixed",
    ,
    drop = FALSE
  ]
  
  if (!is.null(true_beta)) {
    est_ok$true_value <- true_beta[est_ok$term]
    est_ok$bias <- est_ok$estimate - est_ok$true_value
    est_ok$sq_error <- est_ok$bias^2
  } else {
    est_ok$true_value <- NA_real_
    est_ok$bias <- NA_real_
    est_ok$sq_error <- NA_real_
  }
  
  agg_mean <- aggregate(
    estimate ~ method + mechanism + model + term,
    data = est_ok,
    FUN = mean,
    na.rm = TRUE
  )
  
  names(agg_mean)[names(agg_mean) == "estimate"] <- "mean_estimate"
  
  agg_sd <- aggregate(
    estimate ~ method + mechanism + model + term,
    data = est_ok,
    FUN = sd,
    na.rm = TRUE
  )
  
  names(agg_sd)[names(agg_sd) == "estimate"] <- "sd_estimate"
  
  agg_bias <- aggregate(
    bias ~ method + mechanism + model + term,
    data = est_ok,
    FUN = mean,
    na.rm = TRUE
  )
  
  names(agg_bias)[names(agg_bias) == "bias"] <- "mean_bias"
  
  agg_rmse <- aggregate(
    sq_error ~ method + mechanism + model + term,
    data = est_ok,
    FUN = function(x) sqrt(mean(x, na.rm = TRUE))
  )
  
  names(agg_rmse)[names(agg_rmse) == "sq_error"] <- "rmse"
  
  out <- Reduce(
    function(x, y) merge(x, y, by = c("method", "mechanism", "model", "term"), all = TRUE),
    list(agg_mean, agg_sd, agg_bias, agg_rmse)
  )
  
  out[order(out$method, out$mechanism, out$model, out$term), ]
}




extract_bayes_estimates <- function(fit,
                                    fixed_formula,
                                    data,
                                    method = "Bayesian CART-assisted MCMC",
                                    model,
                                    rep,
                                    mechanism,
                                    id_var = "id",
                                    y_var = "y",
                                    ci_probs = c(0.025, 0.975)) {
  
  post <- posterior_draws_after_burnin(fit)
  
  if (is.null(post$beta)) {
    stop("posterior_draws_after_burnin(fit) does not contain post$beta.")
  }
  
  beta_draws <- as.matrix(post$beta)
  
  # Try to recover coefficient names
  X_tmp <- model_matrix_from_formula(fixed_formula, data)
  
  if (ncol(beta_draws) == ncol(X_tmp)) {
    colnames(beta_draws) <- colnames(X_tmp)
  } else if (is.null(colnames(beta_draws))) {
    colnames(beta_draws) <- paste0("beta_", seq_len(ncol(beta_draws)))
  }
  
  q_beta <- apply(beta_draws, 2, stats::quantile, probs = ci_probs, na.rm = TRUE)
  
  out_beta <- data.frame(
    method = method,
    model = model,
    rep = rep,
    mechanism = mechanism,
    term = colnames(beta_draws),
    posterior_mean = colMeans(beta_draws, na.rm = TRUE),
    posterior_median = apply(beta_draws, 2, stats::median, na.rm = TRUE),
    posterior_sd = apply(beta_draws, 2, stats::sd, na.rm = TRUE),
    ci_lower = q_beta[1, ],
    ci_upper = q_beta[2, ],
    component = "fixed",
    row.names = NULL
  )
  
  # Random-effects variance, if available
  sigma_name <- intersect(
    c("sigma2_a", "sigma_a2", "tau2", "sigma2"),
    names(post)
  )
  
  if (length(sigma_name) > 0L) {
    sigma_draws <- as.numeric(post[[sigma_name[1]]])
    
    q_sigma <- stats::quantile(
      sigma_draws,
      probs = ci_probs,
      na.rm = TRUE
    )
    
    out_sigma <- data.frame(
      method = method,
      model = model,
      rep = rep,
      mechanism = mechanism,
      term = "sigma2_a",
      posterior_mean = mean(sigma_draws, na.rm = TRUE),
      posterior_median = median(sigma_draws, na.rm = TRUE),
      posterior_sd = sd(sigma_draws, na.rm = TRUE),
      ci_lower = q_sigma[1],
      ci_upper = q_sigma[2],
      component = "random_effect_variance",
      row.names = NULL
    )
    
    out <- rbind(out_beta, out_sigma)
  } else {
    out <- out_beta
  }
  
  out
}


make_bayes_estimate_table <- function(results,
                                      model_specs,
                                      id_var = "id",
                                      y_var = "y",
                                      output_dir_arg = NULL) {
  
  if (is.null(results$bayes_fits) || length(results$bayes_fits) == 0L) {
    stop("No Bayesian fits found in results$bayes_fits.")
  }
  
  est_list <- list()
  counter <- 1L
  
  for (i in seq_along(results$bayes_fits)) {
    
    sim_item <- results$sim_data[[i]]
    fits_i <- results$bayes_fits[[i]]
    
    rep_id <- sim_item$rep_id
    mech <- sim_item$mechanism
    data_mis <- as.data.frame(sim_item$data_mis)
    
    for (model_name in names(fits_i)) {
      
      message("Bayes estimates | rep=", rep_id,
              " | mechanism=", mech,
              " | model=", model_name)
      
      est_list[[counter]] <- extract_bayes_estimates(
        fit = fits_i[[model_name]],
        fixed_formula = model_specs[[model_name]]$fixed,
        data = data_mis,
        method = "Bayesian CART-assisted MCMC",
        model = model_name,
        rep = rep_id,
        mechanism = mech,
        id_var = id_var,
        y_var = y_var
      )
      
      counter <- counter + 1L
    }
  }
  
  est_tab <- do.call(rbind, est_list)
  
  if (is.null(output_dir_arg)) {
    if (exists("output_dir", envir = .GlobalEnv)) {
      output_dir_arg <- get("output_dir", envir = .GlobalEnv)
    } else {
      output_dir_arg <- "output"
    }
  }
  
  tables_dir <- file.path(output_dir_arg, "tables")
  if (!dir.exists(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE)
  }
  
  utils::write.csv(
    est_tab,
    file.path(tables_dir, "table_bayes_estimates.csv"),
    row.names = FALSE
  )
  
  est_tab
}




read_csv_if_nonempty <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }
  
  info <- file.info(file)
  
  if (is.na(info$size) || info$size == 0L) {
    return(NULL)
  }
  
  out <- tryCatch(
    utils::read.csv(file, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  out
}



read_ml_results <- function(out_dir = file.path("output", "tables")) {
  
  model_file <- file.path(out_dir, "ml_model_results.csv")
  estimate_file <- file.path(out_dir, "ml_estimates.csv")
  selection_file <- file.path(out_dir, "ml_model_selection.csv")
  selection_summary_file <- file.path(out_dir, "ml_model_selection_summary.csv")
  
  list(
    ml = read_csv_if_nonempty(model_file),
    ml_estimates = read_csv_if_nonempty(estimate_file),
    ml_selection = read_csv_if_nonempty(selection_file),
    ml_selection_summary = read_csv_if_nonempty(selection_summary_file)
  )
}



make_estimate_comparison <- function(ml_estimates,
                                     bayes_estimates,
                                     model = "M2_true") {
  true_values <- c(
    "(Intercept)" = 0.20,
    "x2" = -1.50,
    "x3" = 0.80,
    "x4" = -0.30,
    "x5B" = -0.40,
    "x5C" = -0.30,
    "sigma2_a" = 1.00
  )
  
  ml <- subset(ml_estimates, model == model)
  ml$estimate_common <- ml$estimate
  ml$se_common <- ml$std_error
  ml$source <- "ML"
  
  bayes <- subset(bayes_estimates, model == model)
  bayes$estimate_common <- bayes$posterior_mean
  bayes$se_common <- bayes$posterior_sd
  bayes$source <- "Bayes"
  
  keep_cols <- c(
    "source", "method", "model", "rep", "mechanism",
    "term", "estimate_common", "se_common", "component"
  )
  
  common <- rbind(
    ml[, keep_cols],
    bayes[, keep_cols]
  )
  
  common$true_value <- true_values[common$term]
  common$bias <- common$estimate_common - common$true_value
  common$abs_bias <- abs(common$bias)
  
  common[order(common$mechanism, common$term, common$method), ]
}
