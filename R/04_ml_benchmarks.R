# -----------------------------------------------------------------------------
# 04_ml_benchmarks.R
# Complete-case, MICE, missForest, LRT and BIC/log-BF benchmarks
# -----------------------------------------------------------------------------

# measure runtime
with_runtime <- function(expr) {
  t0 <- proc.time()
  value <- force(expr)
  runtime <- as.numeric((proc.time() - t0)[["elapsed"]])
  list(value = value, runtime_seconds = runtime)
}



make_glmer_formula <- function(
  fixed_formula,
  y_var = "y",
  id_var = "id") {
  term_labels <- attr(stats::terms(fixed_formula), "term.labels")
  
  rhs_fixed <- if (length(term_labels) == 0L) {
    "1"
  } else {
    paste(term_labels, collapse = " + ")
  }
  
  f <- stats::as.formula(
    sprintf("%s ~ %s + (1 | %s)", y_var, rhs_fixed, id_var)
  )
  
  environment(f) <- .GlobalEnv
  
  f
}



fit_glmer_safe <- function(
  data, 
  fixed_formula, 
  id_var = "id", 
  y_var = "y",
  nAGQ = 1, 
  maxfun = 2e5
) {
  form <- make_glmer_formula(fixed_formula, id_var = id_var, y_var = y_var)
  tryCatch(
    lme4::glmer(
      formula = form,
      data = data,
      family = binomial(link = "probit"),
      nAGQ = nAGQ,
      control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = maxfun))
    ),
    error = function(e) e
  )
}

extract_ml_fit <- function(
  fit, 
  method, 
  model_name, 
  rep, 
  mechanism = NA_character_, 
  runtime_seconds = NA_real_
) {
  if (inherits(fit, "error")) {
    return(data.frame(
      method = method, model = model_name, rep = rep, mechanism = mechanism,
      converged = FALSE, error = fit$message,
      logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, nobs = NA_integer_,
      runtime_seconds = runtime_seconds
    ))
  }

  data.frame(
    method = method, model = model_name, rep = rep, mechanism = mechanism,
    converged = is.null(fit@optinfo$conv$lme4$messages),
    error = NA_character_,
    logLik = as.numeric(stats::logLik(fit)),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    nobs = stats::nobs(fit),
    runtime_seconds = runtime_seconds
  )
}

fit_complete_case_ml <- function(data, fixed_formula, id_var = "id", y_var = "y", ...) {
  vars <- all.vars(make_glmer_formula(fixed_formula, id_var = id_var, y_var = y_var))
  dat_cc <- stats::na.omit(data[, unique(vars), drop = FALSE])
  with_runtime(fit_glmer_safe(dat_cc, fixed_formula, id_var = id_var, y_var = y_var, ...))
}

fit_mice_ml <- function(
  data,
  fixed_formula,
  id_var = "id",
  y_var = "y",
  m = 5,
  maxit = 5,
  method_default = "cart",
  seed = 123) {
  stopifnot(requireNamespace("mice", quietly = TRUE))
  stopifnot(requireNamespace("lme4", quietly = TRUE))
  
  set.seed(as.integer(seed))
  
  # Extract predictors from fixed-effects formula
  term_labels <- attr(stats::terms(fixed_formula), "term.labels")
  
  # Build glmer formula
  rhs <- if (length(term_labels) == 0L) {
    "1"
  } else {
    paste(term_labels, collapse = " + ")
  }
  
  glmer_formula <- stats::as.formula(
    paste0(y_var, " ~ ", rhs, " + (1 | ", id_var, ")")
  )
  
  data <- as.data.frame(data)
  
  if (!y_var %in% names(data)) {
    stop("fit_mice_ml: y_var '", y_var, "' not found in data. Available columns: ",
         paste(names(data), collapse = ", "))
  }
  
  if (!id_var %in% names(data)) {
    stop("fit_mice_ml: id_var '", id_var, "' not found in data. Available columns: ",
         paste(names(data), collapse = ", "))
  }
  
  data[[id_var]] <- as.factor(data[[id_var]])
  
  # Case 1: Intercept-only model: no covariates to impute
  if (length(term_labels) == 0L) {
    t0 <- proc.time()
    
    fit <- lme4::glmer(
      formula = glmer_formula,
      data = data,
      family = binomial(link = "probit"),
      control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
    
    runtime <- as.numeric((proc.time() - t0)[["elapsed"]])
    
    return(list(
      method = "MICE + ML",
      note = "No fixed-effect covariates in the model; MICE skipped.",
      fits = list(fit),
      pooled = NULL,
      runtime_seconds = runtime,
      logLik = as.numeric(stats::logLik(fit)),
      AIC = stats::AIC(fit),
      BIC = stats::BIC(fit)
    ))
  }
  
  # Variables needed for imputation and model fitting
  vars_needed <- unique(c(y_var, id_var, term_labels))
  vars_needed <- vars_needed[vars_needed %in% names(data)]
  
  dat <- data[, vars_needed, drop = FALSE]
  
  # Ensure outcome and id are available
  if (!all(c(y_var, id_var) %in% names(dat))) {
    stop("Outcome variable or id variable is missing from the data.")
  }
  
  dat[[id_var]] <- as.factor(dat[[id_var]])
  
  # If model covariates have no missing values, skip MICE
  predictor_missing <- any(is.na(dat[, term_labels, drop = FALSE]))
  
  if (!predictor_missing) {
    t0 <- proc.time()
    
    fit <- lme4::glmer(
      formula = glmer_formula,
      data = dat,
      family = binomial(link = "probit"),
      control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
    
    runtime <- as.numeric((proc.time() - t0)[["elapsed"]])
    
    return(list(
      method = "MICE + ML",
      note = "No missing values in model covariates; MICE skipped.",
      fits = list(fit),
      pooled = NULL,
      runtime_seconds = runtime,
      logLik = as.numeric(stats::logLik(fit)),
      AIC = stats::AIC(fit),
      BIC = stats::BIC(fit)
    ))
  }
  
  # MICE setup
  t0 <- proc.time()
  
  ini <- mice::mice(dat, m = m, maxit = 0, printFlag = FALSE)
  pred <- ini$predictorMatrix
  meth <- ini$method
  
  # Do not impute outcome or id
  meth[y_var] <- ""
  meth[id_var] <- ""
  
  # Use default method only for variables with missings
  for (v in names(meth)) {
    if (v %in% c(y_var, id_var)) next
    
    if (any(is.na(dat[[v]]))) {
      meth[v] <- method_default
    } else {
      meth[v] <- ""
    }
  }
  
  # Do not use id as an imputation target/predictor by default
  pred[, id_var] <- 0
  pred[id_var, ] <- 0
  
  # Outcome can be predictor but is not imputed
  pred[y_var, ] <- 0
  
  imp <- mice::mice(
    dat,
    m = m,
    maxit = maxit,
    method = meth,
    predictorMatrix = pred,
    seed = as.integer(seed),
    printFlag = FALSE
  )
  
  completed <- mice::complete(imp, action = "all")
  
  fits <- lapply(completed, function(d) {
    d[[id_var]] <- as.factor(d[[id_var]])
    
    lme4::glmer(
      formula = glmer_formula,
      data = d,
      family = binomial(link = "probit"),
      control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
  })
  
  runtime <- as.numeric((proc.time() - t0)[["elapsed"]])
  
  logLik_values <- sapply(fits, function(f) as.numeric(stats::logLik(f)))
  AIC_values <- sapply(fits, stats::AIC)
  BIC_values <- sapply(fits, stats::BIC)
  
  list(
    method = "MICE + ML",
    note = "MICE imputation used for incomplete model covariates.",
    imp = imp,
    fits = fits,
    runtime_seconds = runtime,
    logLik = median(logLik_values, na.rm = TRUE),
    AIC = median(AIC_values, na.rm = TRUE),
    BIC = median(BIC_values, na.rm = TRUE),
    logLik_values = logLik_values,
    AIC_values = AIC_values,
    BIC_values = BIC_values
  )
}



fit_missforest_ml <- function(
  data,
  fixed_formula,
  id_var = "id",
  y_var = "y",
  ntree = 100,
  maxiter = 10,
  seed = 123,
  ...) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("Package 'missForest' is required for this benchmark.", call. = FALSE)
  }

  set.seed(seed)
  dat <- as.data.frame(data)
  exclude <- c(y_var, id_var)
  xvars <- setdiff(names(dat), exclude)
  xmis <- dat[, xvars, drop = FALSE]
  xmis[] <- lapply(xmis, function(z) if (is.character(z)) factor(z) else z)

  timed <- with_runtime({
    imp <- missForest::missForest(xmis = xmis, maxiter = maxiter, ntree = ntree,
                                  variablewise = TRUE, verbose = FALSE)
    dat_imp <- dat
    dat_imp[, xvars] <- imp$ximp
    fit <- fit_glmer_safe(dat_imp, fixed_formula, id_var = id_var, y_var = y_var, ...)
    list(fit = fit, data_imp = dat_imp, oob_error = imp$OOBerror)
  })

  timed
}

summarise_mice_fits <- function(timed_mice, method, model_name, rep, mechanism = NA_character_) {
  fits <- timed_mice$value$fits
  per_fit <- do.call(rbind, lapply(seq_along(fits), function(i) {
    extract_ml_fit(fits[[i]], method = paste0(method, "_imp", i), model_name = model_name,
                   rep = rep, mechanism = mechanism, runtime_seconds = NA_real_)
  }))

  data.frame(
    method = method,
    model = model_name,
    rep = rep,
    mechanism = mechanism,
    converged = all(per_fit$converged, na.rm = TRUE),
    error = if (all(is.na(per_fit$error))) NA_character_ else paste(na.omit(per_fit$error), collapse = " | "),
    logLik = median(per_fit$logLik, na.rm = TRUE),
    AIC = median(per_fit$AIC, na.rm = TRUE),
    BIC = median(per_fit$BIC, na.rm = TRUE),
    nobs = median(per_fit$nobs, na.rm = TRUE),
    runtime_seconds = timed_mice$runtime_seconds
  )
}

bic_logbf_table <- function(results, reference_model) {
  ref <- results[results$model == reference_model, c("rep", "mechanism", "method", "BIC")]
  names(ref)[names(ref) == "BIC"] <- "BIC_ref"
  out <- merge(results, ref, by = c("rep", "mechanism", "method"), all.x = TRUE)
  out$delta_BIC_vs_ref <- out$BIC - out$BIC_ref
  out$logBF_vs_ref <- -0.5 * out$delta_BIC_vs_ref
  out
}

select_model_by_min_bic <- function(results) {
  split_res <- split(results, list(results$method, results$mechanism, results$rep), drop = TRUE)
  do.call(rbind, lapply(split_res, function(d) {
    d <- d[!is.na(d$BIC), , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d[which.min(d$BIC), c("method", "mechanism", "rep", "model", "BIC")]
  }))
}

pairwise_median_lrt <- function(fits_reduced, fits_full) {
  stopifnot(length(fits_reduced) == length(fits_full))
  pvals <- numeric(length(fits_reduced))
  stat <- numeric(length(fits_reduced))
  df <- numeric(length(fits_reduced))

  for (i in seq_along(fits_reduced)) {
    if (inherits(fits_reduced[[i]], "error") || inherits(fits_full[[i]], "error")) {
      pvals[i] <- NA_real_; stat[i] <- NA_real_; df[i] <- NA_real_
    } else {
      ll0 <- as.numeric(stats::logLik(fits_reduced[[i]]))
      ll1 <- as.numeric(stats::logLik(fits_full[[i]]))
      df0 <- attr(stats::logLik(fits_reduced[[i]]), "df")
      df1 <- attr(stats::logLik(fits_full[[i]]), "df")
      stat[i] <- 2 * (ll1 - ll0)
      df[i] <- df1 - df0
      pvals[i] <- stats::pchisq(stat[i], df = df[i], lower.tail = FALSE)
    }
  }

  list(statistic_median = median(stat, na.rm = TRUE),
       df_median = median(df, na.rm = TRUE),
       p_median = median(pvals, na.rm = TRUE),
       p_values = pvals)
}
