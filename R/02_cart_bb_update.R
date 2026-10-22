# -----------------------------------------------------------------------------
# 02_cart_bb_update.R
# Corrected CART + Bayesian-bootstrap missing covariate update
# -----------------------------------------------------------------------------

# Draw from Bayesian-bootstrap weights over a donor pool.
draw_bb <- function(donor_pool, max_donors = Inf) {
  donor_pool <- donor_pool[!is.na(donor_pool)]
  if (length(donor_pool) == 0) return(NA)
  if (length(donor_pool) == 1) return(donor_pool[[1]])

  if (is.finite(max_donors) && length(donor_pool) > max_donors) {
    donor_pool <- sample(donor_pool, size = max_donors, replace = FALSE)
  }

  w <- rexp(length(donor_pool), rate = 1)
  w <- w / sum(w)
  sample(donor_pool, size = 1, replace = TRUE, prob = w)
}

# This function updates ONLY cells that were originally missing.
# observed cells remain fixed, even after repeated iterations.
cart_bb_update_missing <- function(X_current,
                                   missing_mask,
                                   z = NULL,
                                   a_id = NULL,
                                   y_outcome = NULL,
                                   minbucket_range = 20:60,
                                   maxdepth = 4,
                                   cp = 0.001,
                                   max_donors = 50,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X_current <- as.data.frame(X_current)
  missing_mask <- as.matrix(missing_mask)

  stopifnot(nrow(X_current) == nrow(missing_mask))
  stopifnot(ncol(X_current) == ncol(missing_mask))

  for (k in seq_len(ncol(X_current))) {
    if (!any(missing_mask[, k])) next

    target_name <- names(X_current)[k]
    mis_ind <- which(missing_mask[, k])
    obs_ind <- which(!missing_mask[, k] & !is.na(X_current[[k]]))

    if (length(obs_ind) < 5) next

    y_obs <- X_current[[k]][obs_ind]

    pred_df <- X_current[, -k, drop = FALSE]
    if (!is.null(z)) pred_df$latent_z_current <- as.numeric(z)
    if (!is.null(a_id)) pred_df$random_effect_current <- as.numeric(a_id)
    if (!is.null(y_outcome)) pred_df$outcome_y <- as.numeric(y_outcome)

    df_obs <- data.frame(.target = y_obs, pred_df[obs_ind, , drop = FALSE], check.names = FALSE)
    df_mis <- pred_df[mis_ind, , drop = FALSE]

    # Drop columns that are still entirely missing or have zero variance in the observed training set.
    keep <- vapply(df_obs[, -1, drop = FALSE], function(x) {
      ux <- unique(x[!is.na(x)])
      length(ux) > 1
    }, logical(1))
    if (!any(keep)) next
    df_obs <- data.frame(.target = df_obs$.target, df_obs[, -1, drop = FALSE][, keep, drop = FALSE])
    df_mis <- df_mis[, keep, drop = FALSE]

    is_cat <- is.factor(y_obs) || is.character(y_obs)
    method <- if (is_cat) "class" else "anova"

    if (is_cat) df_obs$.target <- factor(df_obs$.target)

    ctrl <- rpart::rpart.control(
      minbucket = sample(minbucket_range, size = 1),
      maxdepth = maxdepth,
      cp = cp,
      xval = 0
    )

    tree <- rpart::rpart(.target ~ ., data = df_obs, method = method, control = ctrl, model = TRUE)

    if (is_cat) {
      probs <- predict(tree, newdata = df_mis, type = "prob")
      if (is.null(dim(probs))) {
        probs <- matrix(probs, nrow = length(mis_ind), byrow = TRUE)
        colnames(probs) <- levels(df_obs$.target)
      }
      lev <- colnames(probs)
      new_vals <- apply(probs, 1, function(p) sample(lev, size = 1, prob = p / sum(p)))
      if (is.factor(X_current[[k]])) {
        X_current[[k]][mis_ind] <- factor(new_vals, levels = levels(X_current[[k]]))
      } else {
        X_current[[k]][mis_ind] <- new_vals
      }
    } else {
      # Terminal-node proxy: for a regression tree, all observations in one terminal node
      # receive the same prediction. Match missing cases to observed cases with the same
      # terminal-node prediction and draw from those donors.
      pred_obs <- as.numeric(predict(tree, newdata = df_obs[, -1, drop = FALSE], type = "vector"))
      pred_mis <- as.numeric(predict(tree, newdata = df_mis, type = "vector"))

      for (ii in seq_along(mis_ind)) {
        candidates <- which(abs(pred_obs - pred_mis[ii]) < sqrt(.Machine$double.eps))
        if (length(candidates) == 0) {
          # Fallback to nearest terminal prediction.
          candidates <- which.min(abs(pred_obs - pred_mis[ii]))
        }
        donors <- df_obs$.target[candidates]
        X_current[[k]][mis_ind[ii]] <- draw_bb(donors, max_donors = max_donors)
      }
    }
  }

  X_current
}

initialise_missing_values <- function(X, missing_mask = is.na(X), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X0 <- as.data.frame(X)
  missing_mask <- as.matrix(missing_mask)

  for (k in seq_len(ncol(X0))) {
    if (!any(missing_mask[, k])) next
    obs <- X0[[k]][!missing_mask[, k] & !is.na(X0[[k]])]
    if (length(obs) == 0) stop("No observed donors for variable ", names(X0)[k], call. = FALSE)
    if (is.factor(X0[[k]])) {
      X0[[k]][missing_mask[, k]] <- sample(obs, sum(missing_mask[, k]), replace = TRUE)
    } else {
      X0[[k]][missing_mask[, k]] <- sample(obs, sum(missing_mask[, k]), replace = TRUE)
    }
  }
  X0
}
