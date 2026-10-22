# -----------------------------------------------------------------------------
# 09_ml_imputation_estimates.R
# Adds pooled parameter estimates for MICE-CART + ML and missForest + ML.
# Source this file AFTER 04_ml_benchmarks.R and 05_run_simulation.R.
# -----------------------------------------------------------------------------

.null_coalesce <- function(x, y) if (!is.null(x)) x else y

.bind_rows_base <- function(...) {
  xs <- list(...)
  if (length(xs) == 1L && is.list(xs[[1]]) && !is.data.frame(xs[[1]])) xs <- xs[[1]]
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (length(xs) == 0L) return(data.frame())
  all_cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs2 <- lapply(xs, function(x) {
    miss <- setdiff(all_cols, names(x))
    for (cc in miss) x[[cc]] <- NA
    x[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, xs2)
  rownames(out) <- NULL
  out
}

.is_bad_fit <- function(fit) {
  is.null(fit) || inherits(fit, "error") || !inherits(fit, "merMod")
}

.re_var_from_fit <- function(fit, id_var = "id") {
  if (.is_bad_fit(fit)) return(NA_real_)
  vc <- tryCatch(lme4::VarCorr(fit), error = function(e) NULL)
  if (is.null(vc) || length(vc) == 0L) return(NA_real_)
  # Prefer named grouping factor; otherwise first random-effects block.
  block <- if (!is.null(vc[[id_var]])) vc[[id_var]] else vc[[1]]
  sdv <- attr(block, "stddev")
  if (length(sdv) == 0L || !is.finite(sdv[1])) return(NA_real_)
  as.numeric(sdv[1]^2)
}

extract_glmer_estimates <- function(fit,
                                    method,
                                    model_name,
                                    rep,
                                    mechanism = NA_character_,
                                    id_var = "id",
                                    data_version = NA_character_,
                                    impute_mode = NA_character_) {
  if (.is_bad_fit(fit)) {
    msg <- if (inherits(fit, "error")) fit$message else "No valid merMod fit available."
    return(data.frame(
      method = method, model = model_name, rep = rep, mechanism = mechanism,
      data_version = data_version, impute_mode = impute_mode,
      term = NA_character_, estimate = NA_real_, std_error = NA_real_,
      z_value = NA_real_, p_value = NA_real_, component = NA_character_,
      error = msg
    ))
  }

  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(co)) {
    beta <- lme4::fixef(fit)
    se <- sqrt(diag(stats::vcov(fit)))
    co <- cbind(Estimate = beta, `Std. Error` = se,
                `z value` = beta / se,
                `Pr(>|z|)` = 2 * stats::pnorm(-abs(beta / se)))
  }

  fixed <- data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    data_version = data_version,
    impute_mode = impute_mode,
    term = rownames(co),
    estimate = as.numeric(co[, 1]),
    std_error = as.numeric(co[, 2]),
    z_value = if (ncol(co) >= 3) as.numeric(co[, 3]) else as.numeric(co[, 1] / co[, 2]),
    p_value = if (ncol(co) >= 4) as.numeric(co[, 4]) else 2 * stats::pnorm(-abs(co[, 1] / co[, 2])),
    component = "fixed",
    error = NA_character_,
    row.names = NULL
  )

  re <- data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    data_version = data_version,
    impute_mode = impute_mode,
    term = "sigma2_a",
    estimate = .re_var_from_fit(fit, id_var = id_var),
    std_error = NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    component = "random_effect_variance",
    error = NA_character_
  )

  .bind_rows_base(fixed, re)
}

.get_mice_fits <- function(mice_res) {
  if (!is.null(mice_res$value) && !is.null(mice_res$value$fits)) return(mice_res$value$fits)
  if (!is.null(mice_res$fits)) return(mice_res$fits)
  if (!is.null(mice_res$value) && inherits(mice_res$value, "merMod")) return(list(mice_res$value))
  list()
}

pool_mice_glmer_estimates <- function(mice_res,
                                      method,
                                      model_name,
                                      rep,
                                      mechanism = NA_character_,
                                      id_var = "id",
                                      data_version = "mis",
                                      impute_mode = "mice_cart") {
  fits <- .get_mice_fits(mice_res)
  fits <- Filter(function(x) !.is_bad_fit(x), fits)

  if (length(fits) == 0L) {
    return(data.frame(
      method = method, model = model_name, rep = rep, mechanism = mechanism,
      data_version = data_version, impute_mode = impute_mode,
      term = NA_character_, estimate = NA_real_, std_error = NA_real_,
      z_value = NA_real_, p_value = NA_real_, component = NA_character_,
      error = "No valid MICE glmer fits available."
    ))
  }

  betas <- lapply(fits, lme4::fixef)
  all_terms <- unique(unlist(lapply(betas, names), use.names = FALSE))
  B <- do.call(rbind, lapply(betas, function(b) {
    out <- rep(NA_real_, length(all_terms)); names(out) <- all_terms
    out[names(b)] <- as.numeric(b)
    out
  }))

  Vlist <- lapply(fits, function(f) {
    v <- tryCatch(as.matrix(stats::vcov(f)), error = function(e) NULL)
    out <- matrix(NA_real_, nrow = length(all_terms), ncol = length(all_terms),
                  dimnames = list(all_terms, all_terms))
    if (!is.null(v)) {
      common <- intersect(all_terms, rownames(v))
      out[common, common] <- v[common, common, drop = FALSE]
    }
    out
  })

  M <- nrow(B)
  qbar <- colMeans(B, na.rm = TRUE)
  Wdiag <- sapply(seq_along(all_terms), function(j) mean(sapply(Vlist, function(v) v[j, j]), na.rm = TRUE))
  Bdiag <- if (M > 1L) apply(B, 2, stats::var, na.rm = TRUE) else rep(0, length(all_terms))
  Tdiag <- Wdiag + (1 + 1 / M) * Bdiag
  se <- sqrt(Tdiag)
  z <- qbar / se

  fixed <- data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    data_version = data_version,
    impute_mode = impute_mode,
    term = all_terms,
    estimate = as.numeric(qbar),
    std_error = as.numeric(se),
    z_value = as.numeric(z),
    p_value = 2 * stats::pnorm(-abs(z)),
    component = "fixed",
    error = NA_character_,
    row.names = NULL
  )

  s2 <- sapply(fits, .re_var_from_fit, id_var = id_var)
  s2 <- s2[is.finite(s2)]
  re <- data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    data_version = data_version,
    impute_mode = impute_mode,
    term = "sigma2_a",
    estimate = if (length(s2)) mean(s2) else NA_real_,
    std_error = if (length(s2) > 1L) stats::sd(s2) else NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    component = "random_effect_variance",
    error = NA_character_
  )

  .bind_rows_base(fixed, re)
}

.get_sim_mis <- function(sim_item) {
  if (!is.null(sim_item$data_mis)) return(sim_item$data_mis)
  if (!is.null(sim_item$missing)) return(sim_item$missing)
  if (!is.null(sim_item$data)) return(sim_item$data)
  stop("Could not find missing-data component in sim_item.", call. = FALSE)
}

.get_rep <- function(sim_item) .null_coalesce(sim_item$rep_id, .null_coalesce(sim_item$rep, NA_integer_))
.get_mech <- function(sim_item) .null_coalesce(sim_item$mechanism, NA_character_)

# Retained only to document the superseded model-specific imputation workflow.
# The executable workflow uses run_ml_benchmarks_one() defined below.
.run_ml_benchmarks_one_legacy <- function(sim_item,
                                           model_specs,
                                           id_var = "id",
                                           y_var = "y",
                                           mice_m = 20,
                                           mice_maxit = 10,
                                           run_missforest = TRUE,
                                           seed = 123) {
  data_mis <- .get_sim_mis(sim_item)
  rep_id <- .get_rep(sim_item)
  mech <- .get_mech(sim_item)

  data_mis <- as.data.frame(data_mis)
  if (id_var %in% names(data_mis)) data_mis[[id_var]] <- as.factor(data_mis[[id_var]])

  # Common complete-case sample across all candidate models.
  all_model_vars <- unique(unlist(lapply(model_specs, function(x) all.vars(x$fixed)), use.names = FALSE))
  common_cc_vars <- unique(c(y_var, id_var, all_model_vars))
  common_cc_vars <- common_cc_vars[common_cc_vars %in% names(data_mis)]
  cc_rows_common <- stats::complete.cases(data_mis[, common_cc_vars, drop = FALSE])
  data_cc_common <- data_mis[cc_rows_common, , drop = FALSE]
  message("Common complete-case sample size = ", nrow(data_cc_common), " / ", nrow(data_mis))

  common_mice <- tryCatch(
    run_common_mice_cart(
      data_mis, model_specs, id_var = id_var, y_var = y_var,
      m = mice_m, maxit = mice_maxit,
      seed = seed + as.integer(rep_id)
    ),
    error = function(e) e
  )
  mice_fits <- if (inherits(common_mice, "error")) NULL else
    .fit_common_completed_models(
      common_mice$completed, model_specs, id_var = id_var, y_var = y_var
    )

  common_mf <- NULL
  if (isTRUE(run_missforest)) {
    common_mf <- tryCatch(
      run_common_missforest(
        data_mis, model_specs, id_var = id_var, y_var = y_var,
        seed = seed + 1000L + as.integer(rep_id)
      ),
      error = function(e) e
    )
  }

  model_results <- list()
  estimate_results <- list()
  runtime_components <- list()
  counter <- 1L
  ecounter <- 1L
  rcounter <- 1L

  for (model_name in names(model_specs)) {
    message("ML benchmark | model = ", model_name)
    fixed_formula <- model_specs[[model_name]]$fixed

    # Complete cases
    cc <- fit_complete_case_ml(data_cc_common, fixed_formula, id_var = id_var, y_var = y_var)
    model_results[[counter]] <- extract_ml_fit(cc$value, "Complete case + ML", model_name, rep_id, mech, cc$runtime_seconds)
    counter <- counter + 1L
    estimate_results[[ecounter]] <- extract_glmer_estimates(cc$value, "Complete case + ML", model_name, rep_id, mech,
                                                            id_var = id_var, data_version = "mis", impute_mode = "complete_case")
    ecounter <- ecounter + 1L

    # MICE-CART + ML
    mice_res <- fit_mice_ml(data_mis, fixed_formula, id_var = id_var, y_var = y_var,
                            m = mice_m, maxit = mice_maxit, seed = seed + as.integer(rep_id))
    model_results[[counter]] <- summarise_mice_fits(mice_res, "MICE-CART + ML", model_name, rep_id, mech)
    counter <- counter + 1L
    estimate_results[[ecounter]] <- pool_mice_glmer_estimates(mice_res, "MICE-CART + ML", model_name, rep_id, mech,
                                                              id_var = id_var, data_version = "mis", impute_mode = "mice_cart")
    ecounter <- ecounter + 1L

    # missForest + ML
    if (isTRUE(run_missforest)) {
      mf <- tryCatch(
        fit_missforest_ml(data_mis, fixed_formula, id_var = id_var, y_var = y_var,
                          seed = seed + 1000L + as.integer(rep_id)),
        error = function(e) list(value = e, runtime_seconds = NA_real_)
      )
      mf_fit <- if (inherits(mf$value, "error")) mf$value else mf$value$fit
      model_results[[counter]] <- extract_ml_fit(mf_fit, "missForest + ML", model_name, rep_id, mech, mf$runtime_seconds)
      counter <- counter + 1L
      estimate_results[[ecounter]] <- extract_glmer_estimates(mf_fit, "missForest + ML", model_name, rep_id, mech,
                                                              id_var = id_var, data_version = "mis", impute_mode = "missforest")
      ecounter <- ecounter + 1L
    }
  }

  list(
    ml = .bind_rows_base(model_results),
    ml_estimates = .bind_rows_base(estimate_results)
  )
}

run_simulation <- function(n_rep = NULL,
                                    rep_ids = NULL,
                                    mechanisms = c("MCAR", "MAR", "MNAR"),
                                    run_ml = TRUE,
                                    run_bayes = FALSE,
                                    run_missforest = TRUE,
                                    mice_m = 20L,
                                    mice_maxit = 10L,
                                    sim_control_arg = NULL,
                                    model_specs_arg = NULL,
                                    seed = 202605) {
  sim_control_use <- .null_coalesce(sim_control_arg, if (exists("sim_control", inherits = TRUE)) get("sim_control", inherits = TRUE) else NULL)
  model_specs_use <- .null_coalesce(model_specs_arg, if (exists("model_specs", inherits = TRUE)) get("model_specs", inherits = TRUE) else stop("model_specs not found."))

  sim_data <- if (is.null(sim_control_use)) {
    make_simulation_data(n_rep = n_rep, rep_ids = rep_ids, mechanisms = mechanisms)
  } else {
    make_simulation_data(
      n_rep = n_rep, rep_ids = rep_ids, mechanisms = mechanisms,
      sim_control = sim_control_use
    )
  }

  ml_results <- list()
  ml_estimates <- list()
  mice_d2_tests <- list()
  ml_runtime_components <- list()
  bayes_results <- list()
  bayes_fits <- list()

  for (i in seq_along(sim_data)) {
    message("Simulation item ", i, " / ", length(sim_data),
            " | rep=", .get_rep(sim_data[[i]]),
            " | mechanism=", .get_mech(sim_data[[i]]))

    if (isTRUE(run_ml)) {
      tmp_ml <- run_ml_benchmarks_one(
        sim_data[[i]], model_specs_use,
        run_missforest = run_missforest,
        mice_m = mice_m,
        mice_maxit = mice_maxit,
        seed = seed
      )
      ml_results[[i]] <- tmp_ml$ml
      ml_estimates[[i]] <- tmp_ml$ml_estimates
      mice_d2_tests[[i]] <- tmp_ml$mice_d2_tests
      ml_runtime_components[[i]] <- tmp_ml$ml_runtime_components
    }

    if (isTRUE(run_bayes) && exists("run_bayes_one", inherits = TRUE)) {
      br <- run_bayes_one(sim_data[[i]], model_specs_use,
                                   mcmc_control = mcmc_control,
                                   tau2_beta = 1,
                                   seed = seed)
      bayes_results[[i]] <- br$summary
      bayes_fits[[i]] <- br$fits
    }
  }

  out <- list(
    ml = .bind_rows_base(ml_results),
    ml_estimates = .bind_rows_base(ml_estimates),
    mice_d2_tests = .bind_rows_base(mice_d2_tests),
    ml_runtime_components = .bind_rows_base(ml_runtime_components),
    bayes = .bind_rows_base(bayes_results),
    bayes_fits = bayes_fits
  )

  if (exists("output_dir", inherits = TRUE)) {
    out_dir <- get("output_dir", inherits = TRUE)
    dir.create(file.path(out_dir, "rds"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file.path(out_dir, "rds", "simulation_results.rds"))
  }

  out
}

# -----------------------------------------------------------------------------
# Runtime/metadata adjustments for ML model-result rows
# Source this file after 04/05/08. These functions intentionally override
# earlier versions in this file and in 04_ml_benchmarks.R.
# -----------------------------------------------------------------------------

.add_ml_meta <- function(x,
                         data_version = NA_character_,
                         impute_mode = NA_character_) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(x)
  x$data_version <- data_version
  x$impute_mode <- impute_mode
  x
}

.get_runtime_seconds <- function(x) {
  cand <- c(
    if (!is.null(x$runtime_seconds)) x$runtime_seconds else NA_real_,
    if (!is.null(x$value) && !is.null(x$value$runtime_seconds)) x$value$runtime_seconds else NA_real_
  )
  cand <- suppressWarnings(as.numeric(cand))
  cand <- cand[is.finite(cand)]
  if (length(cand) == 0L) NA_real_ else cand[1]
}

summarise_mice_fits <- function(mice_res,
                                method,
                                model_name,
                                rep,
                                mechanism = NA_character_,
                                data_version = "mis",
                                impute_mode = "mice_cart") {
  fits <- .get_mice_fits(mice_res)
  fits <- Filter(function(x) !.is_bad_fit(x), fits)
  runtime <- .get_runtime_seconds(mice_res)

  if (length(fits) == 0L) {
    return(data.frame(
      method = method,
      model = model_name,
      rep = rep,
      mechanism = mechanism,
      data_version = data_version,
      impute_mode = impute_mode,
      converged = FALSE,
      error = "No valid MICE glmer fits available.",
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      nobs = NA_integer_,
      runtime_seconds = runtime
    ))
  }

  per_fit <- lapply(seq_along(fits), function(i) {
    out <- extract_ml_fit(
      fits[[i]],
      method = paste0(method, "_imp", i),
      model_name = model_name,
      rep = rep,
      mechanism = mechanism,
      runtime_seconds = NA_real_
    )
    .add_ml_meta(out, data_version = data_version, impute_mode = impute_mode)
  })

  per_fit <- .bind_rows_base(per_fit)

  if (nrow(per_fit) == 0L) {
    return(data.frame(
      method = method,
      model = model_name,
      rep = rep,
      mechanism = mechanism,
      data_version = data_version,
      impute_mode = impute_mode,
      converged = FALSE,
      error = "MICE fits existed, but no model summaries could be extracted.",
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      nobs = NA_integer_,
      runtime_seconds = runtime
    ))
  }

  ok <- isTRUE(any(per_fit$converged, na.rm = TRUE))

  data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    data_version = data_version,
    impute_mode = impute_mode,
    converged = ok,
    error = if (ok) NA_character_ else paste(unique(na.omit(per_fit$error)), collapse = " | "),
    logLik = stats::median(per_fit$logLik, na.rm = TRUE),
    AIC = stats::median(per_fit$AIC, na.rm = TRUE),
    BIC = stats::median(per_fit$BIC, na.rm = TRUE),
    nobs = suppressWarnings(as.integer(round(stats::median(per_fit$nobs, na.rm = TRUE)))),
    runtime_seconds = runtime
  )
}

run_ml_benchmarks_one <- function(sim_item,
                                           model_specs,
                                           id_var = "id",
                                           y_var = "y",
                                           mice_m = 20,
                                           mice_maxit = 10,
                                           run_missforest = TRUE,
                                           seed = 123) {
  data_mis <- .get_sim_mis(sim_item)
  rep_id <- .get_rep(sim_item)
  mech <- .get_mech(sim_item)

  data_mis <- as.data.frame(data_mis)
  if (id_var %in% names(data_mis)) data_mis[[id_var]] <- as.factor(data_mis[[id_var]])

  all_model_vars <- unique(unlist(lapply(model_specs, function(x) all.vars(x$fixed)), use.names = FALSE))
  common_cc_vars <- unique(c(y_var, id_var, all_model_vars))
  common_cc_vars <- common_cc_vars[common_cc_vars %in% names(data_mis)]
  cc_rows_common <- stats::complete.cases(data_mis[, common_cc_vars, drop = FALSE])
  data_cc_common <- data_mis[cc_rows_common, , drop = FALSE]
  message("Common complete-case sample size = ", nrow(data_cc_common), " / ", nrow(data_mis))

  common_mice <- tryCatch(
    run_common_mice_cart(
      data_mis, model_specs, id_var = id_var, y_var = y_var,
      m = mice_m, maxit = mice_maxit,
      seed = seed + as.integer(rep_id)
    ),
    error = function(e) e
  )
  mice_fits <- if (inherits(common_mice, "error")) NULL else
    .fit_common_completed_models(
      common_mice$completed, model_specs, id_var = id_var, y_var = y_var
    )

  common_mf <- NULL
  if (isTRUE(run_missforest)) {
    common_mf <- tryCatch(
      run_common_missforest(
        data_mis, model_specs, id_var = id_var, y_var = y_var,
        seed = seed + 1000L + as.integer(rep_id)
      ),
      error = function(e) e
    )
  }

  model_results <- list()
  estimate_results <- list()
  runtime_components <- list()
  counter <- 1L
  ecounter <- 1L
  rcounter <- 1L

  for (model_name in names(model_specs)) {
    message("ML benchmark | model = ", model_name)
    fixed_formula <- model_specs[[model_name]]$fixed

    cc <- fit_complete_case_ml(data_cc_common, fixed_formula, id_var = id_var, y_var = y_var)
    model_results[[counter]] <- .add_ml_meta(
      extract_ml_fit(cc$value, "Complete case + ML", model_name, rep_id, mech, cc$runtime_seconds),
      data_version = "mis", impute_mode = "complete_case"
    )
    counter <- counter + 1L
    estimate_results[[ecounter]] <- extract_glmer_estimates(
      cc$value, "Complete case + ML", model_name, rep_id, mech,
      id_var = id_var, data_version = "mis", impute_mode = "complete_case"
    )
    ecounter <- ecounter + 1L

    if (inherits(common_mice, "error")) {
      mice_res <- list(
        fits = list(common_mice),
        runtime_seconds = NA_real_
      )
    } else {
      mice_res <- list(
        fits = mice_fits$fits[[model_name]],
        runtime_seconds = common_mice$runtime_seconds +
          unname(mice_fits$runtime_seconds[[model_name]])
      )
    }
    model_results[[counter]] <- summarise_mice_fits(
      mice_res, "MICE-CART + ML", model_name, rep_id, mech,
      data_version = "mis", impute_mode = "mice_cart"
    )
    counter <- counter + 1L
    estimate_results[[ecounter]] <- pool_mice_glmer_estimates(
      mice_res, "MICE-CART + ML", model_name, rep_id, mech,
      id_var = id_var, data_version = "mis", impute_mode = "mice_cart"
    )
    ecounter <- ecounter + 1L

    runtime_components[[rcounter]] <- data.frame(
      method = "MICE-CART + ML", rep = rep_id, mechanism = mech,
      model = model_name,
      imputation_seconds = if (inherits(common_mice, "error")) NA_real_ else common_mice$runtime_seconds,
      model_fit_seconds = if (is.null(mice_fits)) NA_real_ else unname(mice_fits$runtime_seconds[[model_name]]),
      stringsAsFactors = FALSE
    )
    rcounter <- rcounter + 1L

    if (isTRUE(run_missforest)) {
      mf_fit <- if (inherits(common_mf, "error")) common_mf else common_mf$fits[[model_name]]
      mf_runtime <- if (inherits(common_mf, "error")) NA_real_ else
        common_mf$imputation_runtime_seconds + unname(common_mf$fit_runtime_seconds[[model_name]])
      model_results[[counter]] <- .add_ml_meta(
        extract_ml_fit(mf_fit, "missForest + ML", model_name, rep_id, mech, mf_runtime),
        data_version = "mis", impute_mode = "missforest"
      )
      counter <- counter + 1L
      estimate_results[[ecounter]] <- extract_glmer_estimates(
        mf_fit, "missForest + ML", model_name, rep_id, mech,
        id_var = id_var, data_version = "mis", impute_mode = "missforest"
      )
      ecounter <- ecounter + 1L

      runtime_components[[rcounter]] <- data.frame(
        method = "missForest + ML", rep = rep_id, mechanism = mech,
        model = model_name,
        imputation_seconds = if (inherits(common_mf, "error")) NA_real_ else common_mf$imputation_runtime_seconds,
        model_fit_seconds = if (inherits(common_mf, "error")) NA_real_ else unname(common_mf$fit_runtime_seconds[[model_name]]),
        stringsAsFactors = FALSE
      )
      rcounter <- rcounter + 1L
    }
  }

  d2_tests <- if (inherits(common_mice, "error") || is.null(mice_fits)) {
    data.frame()
  } else {
    pool_mice_d2_tests(mice_fits$fits, rep_id = rep_id, mechanism = mech)
  }

  list(
    ml = .bind_rows_base(model_results),
    ml_estimates = .bind_rows_base(estimate_results),
    mice_d2_tests = d2_tests,
    ml_runtime_components = .bind_rows_base(runtime_components)
  )
}
