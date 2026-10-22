###
# 60_SC6_ml_mice_cart_glmer.R
#
# ML benchmark for the empirical NEPS SC6 illustration.
#
# Workflow:
#   1. Read the prepared Wave-16 analysis data.
#   2. Run one common MICE-CART imputation model on the union of variables in
#      the five candidate models.
#   3. Fit five random-intercept probit GLMMs with lme4::glmer() in every
#      completed dataset.
#   4. Pool fixed effects by Rubin's rules.
#   5. Summarise completed-data logLik, AIC, BIC, random-effect variance.
#   6. Compute nested LRTs per imputation, median-p pooling, and D2 pooling
#      where the required packages support it.
#   7. Compute BIC-based log-Bayes-factor approximations for all pairwise
#      model comparisons on the common imputation scale.
#
# This script is designed to be sourced by 00_master_SC6_empirical_analysis.R
# but can also be run standalone after the processed data have been created.
###

# Settings ----------------------------------------------------------------

options(contrasts = c("contr.treatment", "contr.poly"))

message_ts <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

# Locate project and output folders. When sourced from the master file, these
# objects and environment variables already exist.
if (!exists("PROJECT_DIR", inherits = TRUE)) {
  PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SC6_OUTPUT_DIR <- Sys.getenv(
  "SC6_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "data", "processed")
)

SC6_ML_DIR <- Sys.getenv(
  "SC6_ML_DIR",
  unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "03_ml_benchmark")
)

dir.create(SC6_ML_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("SC6_CONTROL", inherits = TRUE)) {
  SC6_CONTROL <- list(
    outcome = "curr_empl",
    group = "state_id",
    impute_m = 20L,
    impute_maxit = 10L,
    ml_link = "probit",
    ml_nAGQ = 7L,
    seed = 20260628L
  )
}

# Backward-compatible defaults if older master file does not define ml_nAGQ.
if (is.null(SC6_CONTROL$ml_nAGQ)) SC6_CONTROL$ml_nAGQ <- 7L
if (is.null(SC6_CONTROL$impute_m)) SC6_CONTROL$impute_m <- 20L
if (is.null(SC6_CONTROL$impute_maxit)) SC6_CONTROL$impute_maxit <- 10L
if (is.null(SC6_CONTROL$seed)) SC6_CONTROL$seed <- 20260628L
if (is.null(SC6_CONTROL$outcome)) SC6_CONTROL$outcome <- "curr_empl"
if (is.null(SC6_CONTROL$group)) SC6_CONTROL$group <- "state_id"
if (is.null(SC6_CONTROL$ml_link)) SC6_CONTROL$ml_link <- "probit"

set.seed(SC6_CONTROL$seed)

# Package management ------------------------------------------------------

required_packages <- c("dplyr", "tidyr", "purrr", "tibble", "mice", "lme4", "rpart")
optional_packages <- c("mitml", "broom.mixed")

missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0L) {
  stop(
    "Required package(s) missing: ", paste(missing_required, collapse = ", "),
    ". Install them before running the ML benchmark.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(mice)
  library(lme4)
  library(broom.mixed)
  library(mitml)
})

optional_available <- vapply(optional_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(optional_available)) {
  message_ts(
    "Optional package(s) unavailable: ",
    paste(optional_packages[!optional_available], collapse = ", "),
    ". D2 pooling will be skipped if the required support is missing."
  )
}

# Input files -------------------------------------------------------------

prepared_data_rds <- file.path(SC6_OUTPUT_DIR, "SC6_wave16_employment_analysis.rds")
model_variables_rds <- file.path(SC6_OUTPUT_DIR, "SC6_wave16_model_variables.rds")
model_formulas_rds <- file.path(SC6_OUTPUT_DIR, "SC6_wave16_model_formulas.rds")

for (f in c(prepared_data_rds, model_variables_rds, model_formulas_rds)) {
  if (!file.exists(f)) stop("Missing required input file: ", f, call. = FALSE)
}

sc6 <- readRDS(prepared_data_rds)
model_variables <- readRDS(model_variables_rds)
model_formulas <- readRDS(model_formulas_rds)

candidate_models <- paste0("M", 1:5)
if (!all(candidate_models %in% names(model_formulas))) {
  stop("model_formulas must contain M1, M2, M3, M4, and M5.", call. = FALSE)
}
if (!all(candidate_models %in% names(model_variables))) {
  stop("model_variables must contain M1, M2, M3, M4, and M5.", call. = FALSE)
}

outcome_var <- SC6_CONTROL$outcome
group_var <- SC6_CONTROL$group
analysis_vars <- unique(unlist(model_variables[candidate_models], use.names = FALSE))
imputation_vars <- unique(c("ID_t", outcome_var, group_var, analysis_vars))

missing_cols <- setdiff(imputation_vars, names(sc6))
if (length(missing_cols) > 0L) {
  stop("Prepared SC6 data are missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# Data preparation for imputation and GLMM fitting ------------------------

as_binary_numeric <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA_real_, as.numeric(round(x)))
}

coerce_analysis_types <- function(dat, state_levels = NULL, unemp_levels = NULL) {
  dat <- as.data.frame(dat)
  
  if ("ID_t" %in% names(dat)) dat$ID_t <- as.numeric(dat$ID_t)
  
  if (outcome_var %in% names(dat)) {
    dat[[outcome_var]] <- as_binary_numeric(dat[[outcome_var]])
  }
  
  if (group_var %in% names(dat)) {
    if (is.null(state_levels)) state_levels <- levels(factor(dat[[group_var]]))
    dat[[group_var]] <- factor(dat[[group_var]], levels = state_levels)
  }
  
  binary_vars <- intersect(c("female", "migration_background", "father_unemployed", "mother_unemployed"), names(dat))
  for (v in binary_vars) dat[[v]] <- as_binary_numeric(dat[[v]])
  
  continuous_vars <- intersect(c("age_c", "age_c2", "years_education", "ict_wle", "ict_meta", "math_wle", "reading_wle"), names(dat))
  for (v in continuous_vars) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
  
  if ("unemployment_duration_cat" %in% names(dat)) {
    if (is.null(unemp_levels)) {
      unemp_levels <- c("no previous spell", "<1 year", "1 to <2 years", "2 to <5 years", ">=5 years")
      observed_levels <- levels(factor(dat$unemployment_duration_cat))
      if (!all(observed_levels %in% unemp_levels)) unemp_levels <- observed_levels
    }
    dat$unemployment_duration_cat <- factor(dat$unemployment_duration_cat, levels = unemp_levels)
  }
  
  dat
}

state_levels <- levels(factor(sc6[[group_var]]))
unemp_levels <- if ("unemployment_duration_cat" %in% names(sc6)) {
  levels(factor(sc6$unemployment_duration_cat))
} else {
  NULL
}

sc6 <- coerce_analysis_types(sc6, state_levels = state_levels, unemp_levels = unemp_levels)

if (anyNA(sc6[[outcome_var]])) stop("Outcome contains missing values after data preparation.", call. = FALSE)
if (anyNA(sc6[[group_var]])) stop("Grouping variable contains missing values after data preparation.", call. = FALSE)
if (!all(na.omit(sc6[[outcome_var]]) %in% c(0, 1))) {
  stop("Outcome must be coded 0/1 for binomial probit GLMM.", call. = FALSE)
}

imp_data <- sc6[, imputation_vars, drop = FALSE]
imp_data <- coerce_analysis_types(imp_data, state_levels = state_levels, unemp_levels = unemp_levels)

# MICE-CART setup ---------------------------------------------------------

message_ts("Preparing common MICE-CART imputation model.")

method <- mice::make.method(imp_data)
pred <- mice::make.predictorMatrix(imp_data)
where <- is.na(imp_data)

# The ID is neither imputed nor used as a predictor. The outcome and the
# federal-state grouping factor are used as predictors but are not imputed.
non_imputed <- intersect(c("ID_t", outcome_var, group_var), names(imp_data))
method[] <- ""
vars_with_missing <- names(imp_data)[vapply(imp_data, anyNA, logical(1))]
method[setdiff(vars_with_missing, non_imputed)] <- "cart"

pred[, "ID_t"] <- 0L
pred["ID_t", ] <- 0L
for (v in non_imputed) {
  method[v] <- ""
  where[, v] <- FALSE
}

diag(pred) <- 0L

# Do not let deterministic age_c2 be imputed from scratch. It is complete in
# the prepared data by design. If it is not, stop and fix data preparation.
if ("age_c2" %in% names(imp_data) && anyNA(imp_data$age_c2)) {
  stop("age_c2 contains missing values. Recreate age_c and age_c2 in data preparation before imputation.", call. = FALSE)
}

mice_method_tbl <- tibble(variable = names(method), method = unname(method))
mice_pred_tbl <- as.data.frame(pred) %>% tibble::rownames_to_column("target")
write.csv(mice_method_tbl, file.path(SC6_ML_DIR, "mice_cart_methods.csv"), row.names = FALSE)
write.csv(mice_pred_tbl, file.path(SC6_ML_DIR, "mice_cart_predictor_matrix.csv"), row.names = FALSE)

# Pre-imputation audit. CART cannot learn an imputation model for variables
# that are all missing or have only one observed value. In that case, mice logs
# a "constant" event and the completed data may still contain NA values.
pre_mice_audit <- tibble(
  variable = names(imp_data),
  n = nrow(imp_data),
  n_missing = vapply(imp_data, function(x) sum(is.na(x)), integer(1)),
  n_observed = vapply(imp_data, function(x) sum(!is.na(x)), integer(1)),
  n_unique_observed = vapply(imp_data, function(x) length(unique(x[!is.na(x)])), integer(1)),
  observed_values = vapply(imp_data, function(x) {
    vals <- unique(as.character(x[!is.na(x)]))
    vals <- vals[seq_len(min(length(vals), 8L))]
    paste(vals, collapse = " | ")
  }, character(1)),
  method = unname(method[names(imp_data)])
) %>%
  mutate(
    all_missing = .data$n_observed == 0L,
    constant_observed = .data$n_observed > 0L & .data$n_unique_observed < 2L,
    problematic_for_imputation = .data$n_missing > 0L & .data$method != "" &
      (.data$all_missing | .data$constant_observed)
  )
write.csv(pre_mice_audit, file.path(SC6_ML_DIR, "mice_preimputation_variable_audit.csv"), row.names = FALSE)

bad_imp_vars <- pre_mice_audit %>% filter(.data$problematic_for_imputation)
if (nrow(bad_imp_vars) > 0L) {
  stop(
    "At least one incomplete covariate is all-missing or constant among observed cases, ",
    "so CART imputation cannot be learned. Problem variable(s): ",
    paste(bad_imp_vars$variable, collapse = ", "),
    ". See mice_preimputation_variable_audit.csv. Fix the variable construction ",
    "or remove the variable from the prespecified model block before running MICE.",
    call. = FALSE
  )
}

message_ts(
  "Running mice::mice(m = ", SC6_CONTROL$impute_m,
  ", maxit = ", SC6_CONTROL$impute_maxit,
  ", method = CART for incomplete covariates)."
)

imp <- mice::mice(
  data = imp_data,
  m = SC6_CONTROL$impute_m,
  maxit = SC6_CONTROL$impute_maxit,
  method = method,
  predictorMatrix = pred,
  where = where,
  seed = SC6_CONTROL$seed,
  printFlag = TRUE
)

saveRDS(imp, file.path(SC6_ML_DIR, "mice_cart_mids.rds"))

# Basic imputation diagnostics and completed-data summaries ---------------

mice_logged_events <- if (is.null(imp$loggedEvents)) tibble() else as_tibble(imp$loggedEvents)
write.csv(mice_logged_events, file.path(SC6_ML_DIR, "mice_cart_logged_events.csv"), row.names = FALSE)

completed_data <- lapply(seq_len(imp$m), function(i) {
  dat_i <- mice::complete(imp, action = i)
  dat_i <- coerce_analysis_types(dat_i, state_levels = state_levels, unemp_levels = unemp_levels)
  # Ensure no observations are dropped differently across models due to missingness.
  miss_after <- vapply(dat_i[, setdiff(imputation_vars, "ID_t"), drop = FALSE], function(x) sum(is.na(x)), integer(1))
  miss_after <- miss_after[miss_after > 0L]
  if (length(miss_after) > 0L) {
    miss_msg <- paste(paste0(names(miss_after), "=", miss_after), collapse = ", ")
    stop(
      "Completed data set ", i,
      " still contains missing values in analysis variables: ", miss_msg,
      ". Inspect mice_cart_logged_events.csv and mice_preimputation_variable_audit.csv.",
      call. = FALSE
    )
  }
  dat_i
})

saveRDS(completed_data, file.path(SC6_ML_DIR, "mice_cart_completed_data_list.rds"))

completed_miss_check <- bind_rows(lapply(seq_along(completed_data), function(i) {
  completed_data[[i]] %>%
    summarise(across(all_of(setdiff(imputation_vars, "ID_t")), ~ sum(is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
    mutate(.imp = i, .before = 1L)
})) %>% select(.imp, variable, n_missing)
write.csv(completed_miss_check, file.path(SC6_ML_DIR, "completed_data_missing_check.csv"), row.names = FALSE)

# GLMM fitting ------------------------------------------------------------

fit_control <- lme4::glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5),
  calc.derivs = TRUE
)

capture_glmer <- function(formula, data, model_name, imp_index, nAGQ_target) {
  warning_messages <- character(0)
  
  fit_expr <- function(nAGQ_value) {
    withCallingHandlers(
      lme4::glmer(
        formula = formula,
        data = data,
        family = stats::binomial(link = SC6_CONTROL$ml_link),
        nAGQ = nAGQ_value,
        control = fit_control
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }
  
  fit <- tryCatch(fit_expr(nAGQ_target), error = function(e) e)
  actual_nAGQ <- nAGQ_target
  fallback_used <- FALSE
  
  if (inherits(fit, "error") && nAGQ_target != 1L) {
    warning_messages <- c(warning_messages, paste0("nAGQ=", nAGQ_target, " failed: ", conditionMessage(fit)))
    fit <- tryCatch(fit_expr(1L), error = function(e) e)
    actual_nAGQ <- 1L
    fallback_used <- TRUE
  }
  
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      log = tibble(
        model = model_name,
        .imp = imp_index,
        status = "error",
        target_nAGQ = nAGQ_target,
        actual_nAGQ = actual_nAGQ,
        fallback_used = fallback_used,
        singular = NA,
        opt_code = NA_integer_,
        conv_messages = NA_character_,
        warning_messages = paste(unique(warning_messages), collapse = " | "),
        error_message = conditionMessage(fit)
      )
    ))
  }
  
  conv_messages <- fit@optinfo$conv$lme4$messages
  if (is.null(conv_messages)) conv_messages <- NA_character_
  opt_code <- fit@optinfo$conv$opt
  if (is.null(opt_code)) opt_code <- NA_integer_
  
  list(
    fit = fit,
    log = tibble(
      model = model_name,
      .imp = imp_index,
      status = "ok",
      target_nAGQ = nAGQ_target,
      actual_nAGQ = actual_nAGQ,
      fallback_used = fallback_used,
      singular = lme4::isSingular(fit, tol = 1e-4),
      opt_code = as.integer(opt_code),
      conv_messages = paste(conv_messages, collapse = " | "),
      warning_messages = paste(unique(warning_messages), collapse = " | "),
      error_message = NA_character_
    )
  )
}

message_ts("Fitting five random-intercept probit GLMMs in each completed data set.")

fit_objects <- vector("list", length(candidate_models))
names(fit_objects) <- candidate_models
fit_log <- list()

for (mname in candidate_models) {
  message_ts("Fitting ", mname, ": ", deparse(model_formulas[[mname]]))
  fit_objects[[mname]] <- vector("list", length(completed_data))
  
  for (i in seq_along(completed_data)) {
    res <- capture_glmer(
      formula = model_formulas[[mname]],
      data = completed_data[[i]],
      model_name = mname,
      imp_index = i,
      nAGQ_target = SC6_CONTROL$ml_nAGQ
    )
    fit_objects[[mname]][[i]] <- res$fit
    fit_log[[length(fit_log) + 1L]] <- res$log
  }
}

fit_log_tbl <- bind_rows(fit_log)
write.csv(fit_log_tbl, file.path(SC6_ML_DIR, "glmer_fit_log.csv"), row.names = FALSE)

failed_fits <- fit_log_tbl %>% filter(.data$status != "ok")
if (nrow(failed_fits) > 0L) {
  write.csv(failed_fits, file.path(SC6_ML_DIR, "glmer_failed_fits.csv"), row.names = FALSE)
  stop("At least one GLMM fit failed. See glmer_failed_fits.csv.", call. = FALSE)
}

saveRDS(fit_objects, file.path(SC6_ML_DIR, "glmer_fit_objects_by_model.rds"))



# Warning table output
glmer_warning_summary <- fit_log_tbl %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    n_fits = dplyr::n(),
    n_ok = sum(status == "ok"),
    n_fallback = sum(fallback_used),
    n_singular = sum(singular),
    n_conv_messages = sum(!is.na(conv_messages)),
    n_warning_messages = sum(!is.na(warning_messages)),
    .groups = "drop"
  )

readr::write_csv(
  glmer_warning_summary,
  file.path(SC6_ML_DIR, "glmer_warning_summary.csv")
)


# Extraction helpers ------------------------------------------------------

extract_fit_stats <- function(fit, model_name, imp_index) {
  ll <- stats::logLik(fit)
  vc <- as.data.frame(lme4::VarCorr(fit))
  re_var <- vc %>% filter(.data$grp == group_var, .data$var1 == "(Intercept)") %>% pull("vcov")
  if (length(re_var) == 0L) re_var <- NA_real_
  
  tibble(
    model = model_name,
    .imp = imp_index,
    nobs = stats::nobs(fit),
    npar = attr(ll, "df"),
    logLik = as.numeric(ll),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    random_intercept_variance = as.numeric(re_var[1])
  )
}

extract_fixef <- function(fit, model_name, imp_index) {
  beta <- lme4::fixef(fit)
  V <- as.matrix(stats::vcov(fit))
  tibble(
    model = model_name,
    .imp = imp_index,
    term = names(beta),
    estimate = as.numeric(beta),
    variance = as.numeric(diag(V)[names(beta)])
  )
}

fit_stats_long <- bind_rows(lapply(candidate_models, function(mname) {
  bind_rows(lapply(seq_along(fit_objects[[mname]]), function(i) {
    extract_fit_stats(fit_objects[[mname]][[i]], mname, i)
  }))
}))

fixef_long <- bind_rows(lapply(candidate_models, function(mname) {
  bind_rows(lapply(seq_along(fit_objects[[mname]]), function(i) {
    extract_fixef(fit_objects[[mname]][[i]], mname, i)
  }))
}))

write.csv(fit_stats_long, file.path(SC6_ML_DIR, "ml_fit_statistics_by_imputation.csv"), row.names = FALSE)
write.csv(fixef_long, file.path(SC6_ML_DIR, "ml_fixed_effects_by_imputation.csv"), row.names = FALSE)

fit_stats_summary <- fit_stats_long %>%
  group_by(.data$model) %>%
  summarise(
    m = n(),
    nobs = median(.data$nobs),
    npar_median = median(.data$npar),
    logLik_median = median(.data$logLik),
    logLik_mean = mean(.data$logLik),
    logLik_sd = stats::sd(.data$logLik),
    AIC_median = median(.data$AIC),
    AIC_mean = mean(.data$AIC),
    AIC_sd = stats::sd(.data$AIC),
    BIC_median = median(.data$BIC),
    BIC_mean = mean(.data$BIC),
    BIC_sd = stats::sd(.data$BIC),
    re_var_median = median(.data$random_intercept_variance),
    re_var_mean = mean(.data$random_intercept_variance),
    re_var_sd = stats::sd(.data$random_intercept_variance),
    .groups = "drop"
  ) %>%
  arrange(match(.data$model, candidate_models))

write.csv(fit_stats_summary, file.path(SC6_ML_DIR, "ml_fit_statistics_summary.csv"), row.names = FALSE)

# Rubin pooling for fixed effects ----------------------------------------

pool_rubin_one_model <- function(df_model) {
  terms <- unique(df_model$term)
  M <- length(unique(df_model$.imp))
  
  bind_rows(lapply(terms, function(term_i) {
    d <- df_model %>% filter(.data$term == term_i) %>% arrange(.data$.imp)
    if (nrow(d) != M) {
      stop("Term '", term_i, "' was not present in all imputations for model ", unique(df_model$model), call. = FALSE)
    }
    
    q <- d$estimate
    u <- d$variance
    qbar <- mean(q)
    ubar <- mean(u)
    b <- stats::var(q)
    if (is.na(b)) b <- 0
    total_var <- ubar + (1 + 1 / M) * b
    se <- sqrt(total_var)
    riv <- if (isTRUE(all.equal(ubar, 0))) NA_real_ else ((1 + 1 / M) * b) / ubar
    lambda <- if (isTRUE(all.equal(total_var, 0))) NA_real_ else ((1 + 1 / M) * b) / total_var
    
    df_old <- if (is.na(riv) || riv == 0) Inf else (M - 1) * (1 + 1 / riv)^2
    crit <- if (is.infinite(df_old)) stats::qnorm(0.975) else stats::qt(0.975, df = df_old)
    test_stat <- qbar / se
    p_value <- if (is.infinite(df_old)) {
      2 * stats::pnorm(abs(test_stat), lower.tail = FALSE)
    } else {
      2 * stats::pt(abs(test_stat), df = df_old, lower.tail = FALSE)
    }
    
    tibble(
      model = unique(df_model$model),
      term = term_i,
      m = M,
      estimate = qbar,
      std.error = se,
      conf.low = qbar - crit * se,
      conf.high = qbar + crit * se,
      statistic = test_stat,
      df = df_old,
      p.value = p_value,
      ubar = ubar,
      b = b,
      total_var = total_var,
      riv = riv,
      lambda = lambda
    )
  }))
}

coef_pooled <- fixef_long %>%
  group_by(.data$model) %>%
  group_split() %>%
  lapply(pool_rubin_one_model) %>%
  bind_rows() %>%
  arrange(match(.data$model, candidate_models), .data$term)

write.csv(coef_pooled, file.path(SC6_ML_DIR, "ml_coefficients_pooled_rubin.csv"), row.names = FALSE)

# Compact wide coefficient table for inspection. Final journal formatting can
# be handled later in 80_SC6_make_paper_tables.R.
fmt_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

coef_wide <- coef_pooled %>%
  mutate(
    estimate_ci = paste0(fmt_num(.data$estimate), " [", fmt_num(.data$conf.low), ", ", fmt_num(.data$conf.high), "]")
  ) %>%
  select(term, model, estimate_ci) %>%
  tidyr::pivot_wider(names_from = .data$model, values_from = .data$estimate_ci) %>%
  arrange(.data$term)

write.csv(coef_wide, file.path(SC6_ML_DIR, "ml_coefficients_pooled_wide.csv"), row.names = FALSE)

# LRTs: per-imputation, median-p pooling, D2 pooling ----------------------

nested_pairs <- tibble::tribble(
  ~small, ~large, ~comparison,
  "M1", "M2", "Human capital block beyond baseline",
  "M1", "M3", "Labour-market-history block beyond baseline",
  "M1", "M4", "Social-background block beyond baseline",
  "M1", "M5", "Full model beyond baseline",
  "M2", "M5", "Labour-market-history and social-background blocks beyond human capital",
  "M3", "M5", "Human-capital and social-background blocks beyond labour-market history",
  "M4", "M5", "Human-capital and labour-market-history blocks beyond social background"
)

lrt_by_imp <- bind_rows(lapply(seq_len(nrow(nested_pairs)), function(r) {
  small <- nested_pairs$small[r]
  large <- nested_pairs$large[r]
  
  bind_rows(lapply(seq_len(imp$m), function(i) {
    s <- fit_stats_long %>% filter(.data$model == small, .data$.imp == i)
    l <- fit_stats_long %>% filter(.data$model == large, .data$.imp == i)
    df_diff <- l$npar - s$npar
    lr <- 2 * (l$logLik - s$logLik)
    tibble(
      small = small,
      large = large,
      comparison = nested_pairs$comparison[r],
      .imp = i,
      logLik_small = s$logLik,
      logLik_large = l$logLik,
      df_small = s$npar,
      df_large = l$npar,
      df_diff = df_diff,
      LR = lr,
      p_value = stats::pchisq(lr, df = df_diff, lower.tail = FALSE)
    )
  }))
}))

write.csv(lrt_by_imp, file.path(SC6_ML_DIR, "ml_lrt_nested_by_imputation.csv"), row.names = FALSE)

extract_d2_columns <- function(d2_obj) {
  result <- as.data.frame(d2_obj$result)
  if (nrow(result) < 1L) stop("D2 result has no rows.", call. = FALSE)
  nms <- names(result)
  
  find_col <- function(patterns) {
    idx <- unique(unlist(lapply(patterns, function(p) grep(p, nms, ignore.case = TRUE))))
    if (length(idx) == 0L) return(NA_integer_)
    idx[1]
  }
  
  stat_col <- find_col(c("F.value", "statistic", "test", "value"))
  df1_col <- find_col(c("df1", "num"))
  df2_col <- find_col(c("df2", "den"))
  p_col <- find_col(c("p.value", "Pr", "p"))
  riv_col <- find_col(c("riv", "RIV"))
  
  tibble(
    D2_statistic = if (!is.na(stat_col)) as.numeric(result[[stat_col]][1]) else NA_real_,
    D2_df1 = if (!is.na(df1_col)) as.numeric(result[[df1_col]][1]) else NA_real_,
    D2_df2 = if (!is.na(df2_col)) as.numeric(result[[df2_col]][1]) else NA_real_,
    D2_p_value = if (!is.na(p_col)) as.numeric(result[[p_col]][1]) else NA_real_,
    D2_riv = if (!is.na(riv_col)) as.numeric(result[[riv_col]][1]) else NA_real_
  )
}

compute_d2_pair <- function(small, large) {
  if (!requireNamespace("mitml", quietly = TRUE) || !requireNamespace("broom.mixed", quietly = TRUE)) {
    return(tibble(
      small = small,
      large = large,
      D2_statistic = NA_real_,
      D2_df1 = NA_real_,
      D2_df2 = NA_real_,
      D2_p_value = NA_real_,
      D2_riv = NA_real_,
      D2_status = "skipped: requires packages mitml and broom.mixed"
    ))
  }
  
  out <- tryCatch({
    mira_large <- mice::as.mira(fit_objects[[large]])
    mira_small <- mice::as.mira(fit_objects[[small]])
    d2 <- mice::D2(mira_large, mira_small, use = "likelihood")
    cols <- extract_d2_columns(d2)
    bind_cols(tibble(small = small, large = large), cols, tibble(D2_status = "ok"))
  }, error = function(e) {
    tibble(
      small = small,
      large = large,
      D2_statistic = NA_real_,
      D2_df1 = NA_real_,
      D2_df2 = NA_real_,
      D2_p_value = NA_real_,
      D2_riv = NA_real_,
      D2_status = paste0("error: ", conditionMessage(e))
    )
  })
  
  out
}

lrt_summary <- lrt_by_imp %>%
  group_by(.data$small, .data$large, .data$comparison) %>%
  summarise(
    m = n(),
    df_diff_median = median(.data$df_diff),
    LR_median = median(.data$LR),
    LR_mean = mean(.data$LR),
    LR_sd = stats::sd(.data$LR),
    p_median = median(.data$p_value),
    p_min = min(.data$p_value),
    p_max = max(.data$p_value),
    .groups = "drop"
  )

d2_results <- bind_rows(lapply(seq_len(nrow(nested_pairs)), function(r) {
  compute_d2_pair(nested_pairs$small[r], nested_pairs$large[r])
}))

lrt_summary <- lrt_summary %>%
  left_join(d2_results, by = c("small", "large"))

write.csv(lrt_summary, file.path(SC6_ML_DIR, "ml_lrt_nested_pooled.csv"), row.names = FALSE)

# BIC-based log Bayes factors --------------------------------------------

bic_vec <- fit_stats_summary$BIC_median
names(bic_vec) <- fit_stats_summary$model
bic_vec <- bic_vec[candidate_models]

logbf_mat <- outer(
  bic_vec,
  bic_vec,
  FUN = function(bic_i, bic_j) 0.5 * (bic_j - bic_i)
)
rownames(logbf_mat) <- candidate_models
colnames(logbf_mat) <- candidate_models

logbf_tbl <- as.data.frame(logbf_mat) %>%
  tibble::rownames_to_column("model_i")

write.csv(logbf_tbl, file.path(SC6_ML_DIR, "ml_bic_logbf_matrix.csv"), row.names = FALSE)
saveRDS(logbf_mat, file.path(SC6_ML_DIR, "ml_bic_logbf_matrix.rds"))

bic_ranking <- fit_stats_summary %>%
  select(model, BIC_median, AIC_median, logLik_median) %>%
  arrange(.data$BIC_median) %>%
  mutate(rank_BIC = row_number())
write.csv(bic_ranking, file.path(SC6_ML_DIR, "ml_bic_model_ranking.csv"), row.names = FALSE)

# Minimal LaTeX helper outputs -------------------------------------------

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("<", "$<$", x, fixed = TRUE)
  x <- gsub(">", "$>$", x, fixed = TRUE)
  x
}

write_simple_latex_table <- function(df, path, caption, label) {
  df_chr <- df %>% mutate(across(everything(), as.character))
  df_chr[] <- lapply(df_chr, escape_latex)
  align <- paste0("l", paste(rep("r", ncol(df_chr) - 1L), collapse = ""))
  header <- paste(escape_latex(names(df_chr)), collapse = " & ")
  rows <- apply(df_chr, 1, function(z) paste(z, collapse = " & "))
  lines <- c(
    "\\begin{table}[tbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\small",
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, con = path)
}

fit_table_tex <- fit_stats_summary %>%
  transmute(
    Model = .data$model,
    `Median logLik` = fmt_num(.data$logLik_median, 2),
    `Median AIC` = fmt_num(.data$AIC_median, 2),
    `Median BIC` = fmt_num(.data$BIC_median, 2),
    `RE variance` = fmt_num(.data$re_var_median, 3)
  )

lrt_table_tex <- lrt_summary %>%
  transmute(
    Comparison = paste0(.data$small, " vs ", .data$large),
    `Median LR` = fmt_num(.data$LR_median, 2),
    df = fmt_num(.data$df_diff_median, 0),
    `Median p` = fmt_num(.data$p_median, 4),
    `D2 p` = fmt_num(.data$D2_p_value, 4),
    `D2 status` = .data$D2_status
  )

write_simple_latex_table(
  fit_table_tex,
  file.path(SC6_ML_DIR, "table_ml_fit_statistics.tex"),
  caption = "NEPS SC6 ML benchmark: completed-data fit statistics summarised across imputations.",
  label = "tab:sc6_ml_fit_statistics"
)

write_simple_latex_table(
  lrt_table_tex,
  file.path(SC6_ML_DIR, "table_ml_lrt_nested.tex"),
  caption = "NEPS SC6 ML benchmark: nested likelihood-ratio comparisons after MICE-CART.",
  label = "tab:sc6_ml_lrt_nested"
)

# Save compact result object ---------------------------------------------

result_object <- list(
  control = SC6_CONTROL,
  input_files = list(
    prepared_data = prepared_data_rds,
    model_variables = model_variables_rds,
    model_formulas = model_formulas_rds
  ),
  imputation = list(
    mids = imp,
    method = method,
    predictorMatrix = pred,
    loggedEvents = mice_logged_events
  ),
  models = list(
    formulas = model_formulas[candidate_models],
    fit_log = fit_log_tbl,
    fit_stats_long = fit_stats_long,
    fit_stats_summary = fit_stats_summary,
    coefficients_by_imputation = fixef_long,
    coefficients_pooled = coef_pooled,
    lrt_by_imputation = lrt_by_imp,
    lrt_pooled = lrt_summary,
    bic_logbf_matrix = logbf_mat,
    bic_ranking = bic_ranking
  )
)

saveRDS(result_object, file.path(SC6_ML_DIR, "ml_benchmark_results.rds"))

message_ts("ML benchmark finished. Outputs written to: ", SC6_ML_DIR)
message_ts("Primary files:")
message_ts(" - ml_coefficients_pooled_rubin.csv")
message_ts(" - ml_fit_statistics_summary.csv")
message_ts(" - ml_lrt_nested_pooled.csv")
message_ts(" - ml_bic_logbf_matrix.csv")
message_ts(" - ml_benchmark_results.rds")



