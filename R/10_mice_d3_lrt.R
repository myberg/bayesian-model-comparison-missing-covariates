# -----------------------------------------------------------------------------
# 10_mice_d3_lrt.R
# Common ML imputations and D3 pooling for nested model comparisons.
#
# One MICE-CART object and one missForest completion are created per replication
# and mechanism. Every candidate model is then fitted to those same completed
# data, so BIC differences and nested comparisons share the conditioning data.
# -----------------------------------------------------------------------------

.category_from_dummies <- function(b, c, levels = c("A", "B", "C")) {
  b <- suppressWarnings(as.numeric(b))
  c <- suppressWarnings(as.numeric(c))
  out <- rep(NA_character_, length(b))

  out[!is.na(b) & b == 1 & (is.na(c) | c == 0)] <- "B"
  out[!is.na(c) & c == 1 & (is.na(b) | b == 0)] <- "C"
  out[!is.na(b) & b == 0 & !is.na(c) & c == 0] <- "A"

  factor(out, levels = levels)
}

.add_category_factors <- function(data) {
  out <- as.data.frame(data)
  if (all(c("x5B", "x5C") %in% names(out))) {
    out$x5 <- .category_from_dummies(out$x5B, out$x5C)
  }
  if (all(c("x8B", "x8C") %in% names(out))) {
    out$x8 <- .category_from_dummies(out$x8B, out$x8C)
  }
  out
}

.restore_category_dummies <- function(data) {
  out <- as.data.frame(data)
  if ("x5" %in% names(out)) {
    out$x5 <- factor(out$x5, levels = c("A", "B", "C"))
    out$x5B <- as.integer(out$x5 == "B")
    out$x5C <- as.integer(out$x5 == "C")
  }
  if ("x8" %in% names(out)) {
    out$x8 <- factor(out$x8, levels = c("A", "B", "C"))
    out$x8B <- as.integer(out$x8 == "B")
    out$x8C <- as.integer(out$x8 == "C")
  }
  out
}

.common_imputation_frame <- function(data, model_specs_arg,
                                     id_var = "id", y_var = "y") {
  dat <- .add_category_factors(data)
  model_vars <- unique(unlist(lapply(model_specs_arg, function(x) {
    all.vars(x$fixed)
  }), use.names = FALSE))

  imputation_vars <- setdiff(model_vars, c("x5B", "x5C", "x8B", "x8C"))
  if (any(c("x5B", "x5C") %in% model_vars)) imputation_vars <- c(imputation_vars, "x5")
  if (any(c("x8B", "x8C") %in% model_vars)) imputation_vars <- c(imputation_vars, "x8")
  keep <- unique(c(y_var, id_var, imputation_vars))
  keep <- keep[keep %in% names(dat)]

  dat <- dat[, keep, drop = FALSE]
  dat[[id_var]] <- factor(dat[[id_var]])
  dat
}

run_common_mice_cart <- function(data, model_specs_arg,
                                 id_var = "id", y_var = "y",
                                 m = 20L, maxit = 10L, seed = 202605L) {
  stopifnot(requireNamespace("mice", quietly = TRUE))
  dat <- .common_imputation_frame(data, model_specs_arg, id_var, y_var)

  method <- mice::make.method(dat)
  pred <- mice::make.predictorMatrix(dat)
  method[y_var] <- ""
  method[id_var] <- ""
  pred[y_var, ] <- 0
  pred[id_var, ] <- 0
  pred[, id_var] <- 0

  for (v in names(dat)) {
    if (!v %in% c(y_var, id_var)) {
      method[v] <- if (anyNA(dat[[v]])) "cart" else ""
    }
  }

  elapsed <- system.time({
    imp <- mice::mice(
      dat,
      m = as.integer(m),
      maxit = as.integer(maxit),
      method = method,
      predictorMatrix = pred,
      seed = as.integer(seed),
      printFlag = FALSE
    )
    completed <- lapply(mice::complete(imp, action = "all"), function(d) {
      d[[id_var]] <- factor(d[[id_var]])
      .restore_category_dummies(d)
    })
  })[["elapsed"]]

  list(
    imp = imp,
    completed = completed,
    runtime_seconds = as.numeric(elapsed),
    m = as.integer(m),
    maxit = as.integer(maxit)
  )
}

.fit_common_completed_models <- function(completed, model_specs_arg,
                                         id_var = "id", y_var = "y") {
  fits_by_model <- setNames(vector("list", length(model_specs_arg)), names(model_specs_arg))
  runtime_by_model <- setNames(numeric(length(model_specs_arg)), names(model_specs_arg))

  for (model_name in names(model_specs_arg)) {
    fixed_formula <- model_specs_arg[[model_name]]$fixed
    elapsed <- system.time({
      fits <- lapply(completed, function(d) {
        d[[id_var]] <- factor(d[[id_var]])
        fit_glmer_safe(d, fixed_formula, id_var = id_var, y_var = y_var)
      })
    })[["elapsed"]]
    fits_by_model[[model_name]] <- fits
    runtime_by_model[[model_name]] <- as.numeric(elapsed)
  }

  list(fits = fits_by_model, runtime_seconds = runtime_by_model)
}

run_common_missforest <- function(data, model_specs_arg,
                                  id_var = "id", y_var = "y",
                                  ntree = 100L, maxiter = 10L,
                                  seed = 203605L) {
  stopifnot(requireNamespace("missForest", quietly = TRUE))
  dat <- .common_imputation_frame(data, model_specs_arg, id_var, y_var)
  set.seed(as.integer(seed))

  # The fully observed outcome remains in xmis and can predict incomplete X.
  xmis <- dat[, setdiff(names(dat), id_var), drop = FALSE]
  elapsed_imp <- system.time({
    imp <- missForest::missForest(
      xmis = xmis,
      ntree = as.integer(ntree),
      maxiter = as.integer(maxiter),
      variablewise = TRUE,
      verbose = FALSE
    )
  })[["elapsed"]]

  completed <- dat
  completed[, names(imp$ximp)] <- imp$ximp
  completed[[id_var]] <- factor(completed[[id_var]])
  completed <- .restore_category_dummies(completed)
  fits <- .fit_common_completed_models(
    list(completed), model_specs_arg, id_var = id_var, y_var = y_var
  )

  list(
    completed = completed,
    fits = lapply(fits$fits, function(z) z[[1]]),
    imputation_runtime_seconds = as.numeric(elapsed_imp),
    fit_runtime_seconds = fits$runtime_seconds,
    oob_error = imp$OOBerror
  )
}

default_nested_model_pairs <- function() {
  data.frame(
    null_model = c("M1_intercept", "M3_reduced", "M2_true", "M3_reduced", "M4_noise"),
    alt_model = c("M2_true", "M2_true", "M5_full", "M6_alt", "M5_full"),
    comparison = c("M1 vs M2", "M3 vs M2", "M2 vs M5", "M3 vs M6", "M4 vs M5"),
    stringsAsFactors = FALSE
  )
}

pool_mice_d3_tests <- function(fits_by_model, rep_id, mechanism,
                               nested_pairs = default_nested_model_pairs()) {
  stopifnot(requireNamespace("mitml", quietly = TRUE))
  out <- vector("list", nrow(nested_pairs))

  for (i in seq_len(nrow(nested_pairs))) {
    pair <- nested_pairs[i, ]
    full <- fits_by_model[[pair$alt_model]]
    restricted <- fits_by_model[[pair$null_model]]
    valid <- vapply(full, inherits, logical(1), what = "merMod") &
      vapply(restricted, inherits, logical(1), what = "merMod")

    if (!all(valid)) {
      out[[i]] <- data.frame(
        method = "MICE-CART + ML", rep = rep_id, mechanism = mechanism,
        comparison = pair$comparison, null_model = pair$null_model,
        alt_model = pair$alt_model, pooling = "D3", F_value = NA_real_,
        df1 = NA_real_, df2 = NA_real_, p_value = NA_real_, riv = NA_real_,
        n_imp = sum(valid),
        error = paste0(
          "D3 not evaluated: matched model fits were available for ",
          sum(valid), " of ", length(valid), " imputations."
        ),
        stringsAsFactors = FALSE
      )
      next
    }

    test <- tryCatch(
      mitml::testModels(full[valid], restricted[valid], method = "D3"),
      error = function(e) e
    )
    if (inherits(test, "error")) {
      out[[i]] <- data.frame(
        method = "MICE-CART + ML", rep = rep_id, mechanism = mechanism,
        comparison = pair$comparison, null_model = pair$null_model,
        alt_model = pair$alt_model, pooling = "D3", F_value = NA_real_,
        df1 = NA_real_, df2 = NA_real_, p_value = NA_real_, riv = NA_real_,
        n_imp = sum(valid), error = conditionMessage(test),
        stringsAsFactors = FALSE
      )
    } else {
      stat <- test$test[1, ]
      out[[i]] <- data.frame(
        method = "MICE-CART + ML", rep = rep_id, mechanism = mechanism,
        comparison = pair$comparison, null_model = pair$null_model,
        alt_model = pair$alt_model, pooling = "D3",
        F_value = unname(stat[["F.value"]]),
        df1 = unname(stat[["df1"]]),
        df2 = unname(stat[["df2"]]),
        p_value = unname(stat[["P(>F)"]]),
        riv = unname(stat[["RIV"]]),
        n_imp = sum(valid), error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, out)
}
