# -----------------------------------------------------------------------------
# 08_ml_before_deletion.R
# Adds the full-data / before-deletion maximum-likelihood benchmark.
#
# This is deliberately implemented as a non-invasive add-on: it reads the
# exported full-data CSV files from data/simulated/simulation_index.csv and fits
# the same glmer/probit random-effects models as the other ML benchmarks.
# -----------------------------------------------------------------------------

.resolve_index_path <- function(p, index_file) {
  if (is.na(p) || !nzchar(p)) return(NA_character_)
  if (file.exists(p)) return(p)
  local_path <- file.path(dirname(index_file), basename(p))
  if (file.exists(local_path)) return(local_path)
  p
}

.bind_rows_base <- function(...) {
  xs <- list(...)
  
  # Falls genau eine Liste übergeben wurde, diese entpacken
  if (length(xs) == 1L && is.list(xs[[1]]) && !is.data.frame(xs[[1]])) {
    xs <- xs[[1]]
  }
  
  if (length(xs) == 0L) {
    return(data.frame())
  }
  
  # Nur echte data.frames mit mindestens einer Zeile behalten
  xs <- Filter(
    function(x) {
      is.data.frame(x) && nrow(x) > 0L
    },
    xs
  )
  
  if (length(xs) == 0L) {
    return(data.frame())
  }
  
  all_cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  
  xs2 <- lapply(xs, function(x) {
    missing_cols <- setdiff(all_cols, names(x))
    
    for (cc in missing_cols) {
      x[[cc]] <- NA
    }
    
    x[, all_cols, drop = FALSE]
  })
  
  out <- do.call(rbind, xs2)
  rownames(out) <- NULL
  out
}



.extract_ml_fit_fallback <- function(fit,
                                     method,
                                     model_name,
                                     rep,
                                     mechanism = NA_character_,
                                     runtime_seconds = NA_real_,
                                     data_version = "full",
                                     impute_mode = "none") {
  if (inherits(fit, "error")) {
    return(data.frame(
      method = method,
      model = model_name,
      rep = rep,
      mechanism = mechanism,
      data_version = data_version,
      impute_mode = impute_mode,
      converged = FALSE,
      error = conditionMessage(fit),
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      nobs = NA_integer_,
      runtime_seconds = runtime_seconds,
      stringsAsFactors = FALSE
    ))
  }

  conv_msg <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)

  data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    data_version = data_version,
    impute_mode = impute_mode,
    converged = is.null(conv_msg),
    error = NA_character_,
    logLik = as.numeric(stats::logLik(fit)),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    nobs = stats::nobs(fit),
    runtime_seconds = runtime_seconds,
    stringsAsFactors = FALSE
  )
}

.extract_glmer_estimates_fallback <- function(fit,
                                              method,
                                              model_name,
                                              rep,
                                              mechanism = NA_character_,
                                              data_version = "full",
                                              impute_mode = "none") {
  if (inherits(fit, "error")) {
    return(data.frame(
      method = method,
      model = model_name,
      rep = rep,
      mechanism = mechanism,
      data_version = data_version,
      impute_mode = impute_mode,
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      component = NA_character_,
      error = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)

  fixed <- if (!is.null(co) && nrow(co) > 0L) {
    z_col <- grep("z", colnames(co), ignore.case = TRUE, value = TRUE)[1L]
    p_col <- grep("Pr", colnames(co), ignore.case = TRUE, value = TRUE)[1L]
    data.frame(
      method = method,
      model = model_name,
      rep = rep,
      mechanism = mechanism,
      data_version = data_version,
      impute_mode = impute_mode,
      term = rownames(co),
      estimate = co[, "Estimate"],
      std_error = co[, "Std. Error"],
      z_value = if (!is.na(z_col)) co[, z_col] else NA_real_,
      p_value = if (!is.na(p_col)) co[, p_col] else NA_real_,
      component = "fixed",
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    beta <- lme4::fixef(fit)
    se <- sqrt(diag(stats::vcov(fit)))
    z <- beta / se
    data.frame(
      method = method,
      model = model_name,
      rep = rep,
      mechanism = mechanism,
      data_version = data_version,
      impute_mode = impute_mode,
      term = names(beta),
      estimate = as.numeric(beta),
      std_error = as.numeric(se),
      z_value = as.numeric(z),
      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      component = "fixed",
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  vc <- tryCatch(as.data.frame(lme4::VarCorr(fit)), error = function(e) NULL)
  sigma2_a <- if (!is.null(vc) && "vcov" %in% names(vc) && nrow(vc) > 0L) vc$vcov[1L] else NA_real_

  ran <- data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    data_version = data_version,
    impute_mode = impute_mode,
    term = "sigma2_a",
    estimate = sigma2_a,
    std_error = NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    component = "random_effect_variance",
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  rbind(fixed, ran)
}

.extract_ml_fit_with_meta <- function(fit,
                                      method,
                                      model_name,
                                      rep,
                                      mechanism,
                                      runtime_seconds,
                                      data_version = "full",
                                      impute_mode = "none") {
  if (exists("extract_ml_fit", mode = "function")) {
    out <- extract_ml_fit(
      fit = fit,
      method = method,
      model_name = model_name,
      rep = rep,
      mechanism = mechanism,
      runtime_seconds = runtime_seconds
    )
  } else {
    out <- .extract_ml_fit_fallback(
      fit = fit,
      method = method,
      model_name = model_name,
      rep = rep,
      mechanism = mechanism,
      runtime_seconds = runtime_seconds,
      data_version = data_version,
      impute_mode = impute_mode
    )
  }

  out$data_version <- data_version
  out$impute_mode <- impute_mode
  out
}

.extract_glmer_estimates_with_meta <- function(fit,
                                               method,
                                               model_name,
                                               rep,
                                               mechanism,
                                               data_version = "full",
                                               impute_mode = "none") {
  out <- NULL

  if (exists("extract_glmer_estimates", mode = "function")) {
    out <- tryCatch(
      extract_glmer_estimates(
        fit = fit,
        method = method,
        model_name = model_name,
        rep = rep,
        mechanism = mechanism
      ),
      error = function(e) NULL
    )
  }

  if (is.null(out)) {
    out <- .extract_glmer_estimates_fallback(
      fit = fit,
      method = method,
      model_name = model_name,
      rep = rep,
      mechanism = mechanism,
      data_version = data_version,
      impute_mode = impute_mode
    )
  }

  out$data_version <- data_version
  out$impute_mode <- impute_mode
  out
}

# Main user-facing function.
run_ml_before_deletion_from_index <- function(index_file = file.path("data", "simulated", "simulation_index.csv"),
                                             model_specs_arg = NULL,
                                             out_dir = file.path("output", "tables"),
                                             method_label = "Before deletion + ML",
                                             id_var = "id",
                                             y_var = "y",
                                             models = NULL,
                                             reps = NULL,
                                             mechanisms_filter = NULL,
                                             nAGQ = 1,
                                             maxfun = 2e5,
                                             write_csv = TRUE) {
  if (!file.exists(index_file)) {
    stop("simulation_index.csv not found: ", index_file,
         ". Run export_simulation_data() first.")
  }
  if (!exists("fit_glmer_safe", mode = "function")) {
    stop("fit_glmer_safe() not found. Source R/04_ml_benchmarks.R first.")
  }

  if (is.null(model_specs_arg)) {
    if (exists("model_specs", envir = .GlobalEnv)) {
      model_specs_arg <- get("model_specs", envir = .GlobalEnv)
    } else {
      stop("No model_specs supplied and no global object 'model_specs' found.")
    }
  }

  if (!is.null(models)) {
    miss_models <- setdiff(models, names(model_specs_arg))
    if (length(miss_models) > 0L) {
      stop("Unknown model(s): ", paste(miss_models, collapse = ", "))
    }
    model_specs_arg <- model_specs_arg[models]
  }

  idx <- utils::read.csv(index_file, stringsAsFactors = FALSE)
  if (!is.null(reps)) idx <- idx[idx$rep_id %in% as.integer(reps), , drop = FALSE]
  if (!is.null(mechanisms_filter)) idx <- idx[idx$mechanism %in% mechanisms_filter, , drop = FALSE]
  if (nrow(idx) == 0L) stop("No index rows remain after filtering.")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  model_results <- list()
  estimate_results <- list()
  counter <- 1L

  for (ii in seq_len(nrow(idx))) {
    rep_id <- idx$rep_id[ii]
    mech <- idx$mechanism[ii]
    full_file <- .resolve_index_path(idx$file_full[ii], index_file)

    message("ML before deletion | rep=", rep_id, " | ", mech,
            " | file=", basename(full_file))

    dat_full <- utils::read.csv(full_file, stringsAsFactors = FALSE)
    if (id_var %in% names(dat_full)) dat_full[[id_var]] <- as.factor(dat_full[[id_var]])

    for (model_name in names(model_specs_arg)) {
      fixed_formula <- model_specs_arg[[model_name]]$fixed

      timed <- with_runtime({
        fit_glmer_safe(
          data = dat_full,
          fixed_formula = fixed_formula,
          id_var = id_var,
          y_var = y_var,
          nAGQ = nAGQ,
          maxfun = maxfun
        )
      })

      model_results[[counter]] <- .extract_ml_fit_with_meta(
        fit = timed$value,
        method = method_label,
        model_name = model_name,
        rep = rep_id,
        mechanism = mech,
        runtime_seconds = timed$runtime_seconds,
        data_version = "full",
        impute_mode = "none"
      )

      estimate_results[[counter]] <- .extract_glmer_estimates_with_meta(
        fit = timed$value,
        method = method_label,
        model_name = model_name,
        rep = rep_id,
        mechanism = mech,
        data_version = "full",
        impute_mode = "none"
      )

      counter <- counter + 1L
    }
  }

  out <- list(
    ml = .bind_rows_base(model_results),
    ml_estimates = .bind_rows_base(estimate_results)
  )

  if (isTRUE(write_csv)) {
    utils::write.csv(out$ml,
                     file.path(out_dir, "ml_before_deletion_model_results.csv"),
                     row.names = FALSE)
    utils::write.csv(out$ml_estimates,
                     file.path(out_dir, "ml_before_deletion_estimates.csv"),
                     row.names = FALSE)
  }

  out
}

append_ml_before_deletion <- function(ml_res, before_deletion_res) {
  if (is.null(ml_res)) ml_res <- list(ml = NULL, ml_estimates = NULL)
  ml_res$ml <- .bind_rows_base(ml_res$ml, before_deletion_res$ml)
  ml_res$ml_estimates <- .bind_rows_base(ml_res$ml_estimates, before_deletion_res$ml_estimates)
  ml_res
}
