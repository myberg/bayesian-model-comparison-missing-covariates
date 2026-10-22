# -----------------------------------------------------------------------------
# 03_bayes_sampler.R
# Bayesian random-effects probit sampler with repeated CART updates
# -----------------------------------------------------------------------------

make_fixed_formula <- function(formula_with_y) {
  trm <- stats::terms(formula_with_y)
  term_labels <- attr(trm, "term.labels")
  
  if (length(term_labels) == 0L) {
    return(stats::as.formula("~ 1"))
  }
  
  out <- stats::reformulate(term_labels)
  environment(out) <- environment(formula_with_y)
  out
}

model_matrix_from_formula <- function(fixed_formula, data) {
  if (length(fixed_formula) == 3L) {
    fixed_formula <- make_fixed_formula(fixed_formula)
  }
  
  stats::model.matrix(
    fixed_formula,
    data = data,
    na.action = stats::na.pass
  )
}

with_runtime <- function(expr) {
  t0 <- proc.time()
  value <- force(expr)
  runtime <- as.numeric((proc.time() - t0)[["elapsed"]])
  list(value = value, runtime_seconds = runtime)
}

fit_bayes_re_probit_cart <- function(
  data,
  fixed_formula,
  id_var = "id",
  y_var = "y",
  n_iter = 4000,
  burnin = 1000,
  beta_prior = NULL,
  c0 = 2,
  d0 = 1,
  cart_control = list(
    minbucket_range = 20:60,
    maxdepth = 4,
    cp = 0.001,
    max_donors = 50
  ),
  seed = NULL,
  verbose_every = 1000) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(y_var %in% names(data), id_var %in% names(data))
  dat <- as.data.frame(data)
  dat[[id_var]] <- factor(dat[[id_var]])
  y <- as.integer(dat[[y_var]])
  id <- as.integer(dat[[id_var]])
  J <- length(levels(dat[[id_var]]))
  N <- nrow(dat)
  group_indices <- split(seq_len(N), id)

  # Determine covariates used for missingness updates: all variables excluding outcome and id.
  X_cov <- dat[, setdiff(names(dat), c(y_var, id_var, "z_true", "a_true")), drop = FALSE]
  missing_mask <- is.na(X_cov)
  X_current <- initialise_missing_values(X_cov, missing_mask, seed = seed)

  # Initial design matrix determines p and column names.
  design_data <- cbind(dat[, c(y_var, id_var), drop = FALSE], X_current)
  X_design <- model_matrix_from_formula(fixed_formula, design_data)
  p <- ncol(X_design)

  if (is.null(beta_prior)) beta_prior <- make_beta_prior(p = p, tau2_beta = 1)
  if (length(beta_prior$mu) != p) stop("beta_prior dimension does not match design matrix", call. = FALSE)

  beta <- rep(0, p)
  sigma2_a <- 1
  a <- rep(0, J)
  z <- ifelse(y == 1, 0.5, -0.5)

  beta_store <- matrix(NA_real_, nrow = n_iter, ncol = p)
  colnames(beta_store) <- colnames(X_design)
  sigma2_a_store <- numeric(n_iter)
  a_store <- matrix(NA_real_, nrow = n_iter, ncol = J)
  z_store <- matrix(NA_real_, nrow = n_iter, ncol = N)
  m_beta_store <- matrix(NA_real_, nrow = n_iter, ncol = p)
  V_beta_store <- vector("list", n_iter)
  log_lik_store <- numeric(n_iter)
  X_store <- vector("list", n_iter)

  Sigma_inv <- beta_prior$Sigma_inv
  mu_prior <- beta_prior$mu

  for (iter in seq_len(n_iter)) {
    if (!is.null(verbose_every) && verbose_every > 0 && iter %% verbose_every == 0) {
      message("Bayes CART iteration ", iter, " / ", n_iter)
    }

    # 1. CART-assisted update of the originally missing covariate cells.
    X_current <- cart_bb_update_missing(
      X_current = X_current,
      missing_mask = missing_mask,
      z = z,
      a_id = a[id],
      y_outcome = y,
      minbucket_range = cart_control$minbucket_range,
      maxdepth = cart_control$maxdepth,
      cp = cart_control$cp,
      max_donors = cart_control$max_donors
    )

    design_data <- cbind(dat[, c(y_var, id_var), drop = FALSE], X_current)
    X_design <- model_matrix_from_formula(fixed_formula, design_data)

    # 2. Latent Gaussian responses.
    mu <- as.numeric(X_design %*% beta + a[id])
    z[y == 1] <- truncnorm::rtruncnorm(sum(y == 1), a = 0, b = Inf, mean = mu[y == 1], sd = 1)
    z[y == 0] <- truncnorm::rtruncnorm(sum(y == 0), a = -Inf, b = 0, mean = mu[y == 0], sd = 1)

    # 3. Regression coefficients.
    V_beta <- solve(Sigma_inv + crossprod(X_design))
    m_beta <- V_beta %*% (Sigma_inv %*% mu_prior + crossprod(X_design, z - a[id]))
    beta <- as.numeric(mvtnorm::rmvnorm(1, mean = as.numeric(m_beta), sigma = V_beta))

    # 4. Random intercepts.
    for (j in seq_len(J)) {
      idx <- group_indices[[j]]
      V_a <- 1 / (length(idx) + 1 / sigma2_a)
      m_a <- V_a * sum(z[idx] - as.numeric(X_design[idx, , drop = FALSE] %*% beta))
      a[j] <- rnorm(1, mean = m_a, sd = sqrt(V_a))
    }

    # 5. Random-effect variance.
    shape <- c0 + J / 2
    rate <- d0 + 0.5 * sum(a^2)
    sigma2_a <- 1 / rgamma(1, shape = shape, rate = rate)

    # Store.
    beta_store[iter, ] <- beta
    sigma2_a_store[iter] <- sigma2_a
    a_store[iter, ] <- a
    z_store[iter, ] <- z
    m_beta_store[iter, ] <- as.numeric(m_beta)
    V_beta_store[[iter]] <- V_beta
    X_store[[iter]] <- X_current

    eta <- as.numeric(X_design %*% beta + a[id])
    log_lik_store[iter] <- sum(stats::dbinom(y, size = 1, prob = stats::pnorm(eta), log = TRUE))
  }

  draws <- list(
    beta = beta_store,
    sigma2_a = sigma2_a_store,
    a = a_store,
    z = z_store,
    m_beta = m_beta_store,
    V_beta = V_beta_store,
    log_lik = log_lik_store,
    X_mis = X_store,
    missing_mask = missing_mask
  )

  structure(list(
    draws = draws,
    fixed_formula = fixed_formula,
    id_var = id_var,
    y_var = y_var,
    beta_prior = beta_prior,
    c0 = c0,
    d0 = d0,
    n_iter = n_iter,
    burnin = burnin,
    cart_control = cart_control
  ), class = "smap_bayes_cart_fit")
}

posterior_draws_after_burnin <- function(fit) {
  keep <- seq.int(fit$burnin + 1, fit$n_iter)
  list(
    beta = fit$draws$beta[keep, , drop = FALSE],
    sigma2_a = fit$draws$sigma2_a[keep],
    a = fit$draws$a[keep, , drop = FALSE],
    z = fit$draws$z[keep, , drop = FALSE],
    log_lik = fit$draws$log_lik[keep]
  )
}

summarise_bayes_fit <- function(fit) {
  d <- posterior_draws_after_burnin(fit)
  beta_summary <- t(apply(d$beta, 2, function(x) {
    c(mean = mean(x), median = median(x), sd = sd(x),
      q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)))
  }))

  data.frame(
    parameter = rownames(beta_summary),
    beta_summary,
    row.names = NULL
  )
}

summarise_mcmc_diagnostics <- function(fit) {
  d <- posterior_draws_after_burnin(fit)
  beta_ess <- coda::effectiveSize(coda::mcmc(d$beta))
  sig_ess <- coda::effectiveSize(coda::mcmc(d$sigma2_a))

  data.frame(
    parameter = c(names(beta_ess), "sigma2_a"),
    ess = c(as.numeric(beta_ess), as.numeric(sig_ess)),
    stringsAsFactors = FALSE
  )
}

get_chib_fixed_points <- function(fit) {
  d <- posterior_draws_after_burnin(fit)
  beta_mean <- colMeans(d$beta)
  beta_median <- apply(d$beta, 2, median)
  idx_hp <- which.max(d$log_lik)

  list(
    posterior_mean = list(
      beta = beta_mean, a = colMeans(d$a), sigma2_a = mean(d$sigma2_a)
    ),
    posterior_median = list(
      beta = beta_median, a = apply(d$a, 2, median),
      sigma2_a = median(d$sigma2_a)
    ),
    high_density_draw = list(
      beta = d$beta[idx_hp, ], a = d$a[idx_hp, ],
      sigma2_a = d$sigma2_a[idx_hp]
    )
  )
}

# Select one common completion from a reference CART-BB run. Candidate models
# must all be refitted on the returned data set for conditional Bayes factors.
select_medoid_completed_data <- function(fit, data, medoid_thin = 10L) {
  keep <- seq.int(fit$burnin + 1L, fit$n_iter, by = as.integer(medoid_thin))
  mask <- fit$draws$missing_mask
  covariate_names <- colnames(mask)
  if (is.null(covariate_names)) {
    covariate_names <- setdiff(names(data), c(fit$y_var, fit$id_var, "z_true", "a_true"))
  }

  if (!any(mask)) {
    return(list(data = data, draw_id = keep[1L], total_distance = 0,
                n_draws = length(keep), n_missing_cells = 0L))
  }

  encode_column <- function(x) {
    if (is.factor(x)) return(as.numeric(x))
    if (is.logical(x)) return(as.numeric(x))
    as.numeric(x)
  }
  extract_missing <- function(x) {
    unlist(lapply(seq_along(covariate_names), function(j) {
      encode_column(x[[covariate_names[j]]])[mask[, j]]
    }), use.names = FALSE)
  }

  draw_matrix <- t(vapply(
    fit$draws$X_mis[keep], extract_missing,
    numeric(sum(mask))
  ))

  original_covariates <- data[, covariate_names, drop = FALSE]
  scale_by_column <- vapply(seq_along(covariate_names), function(j) {
    observed <- encode_column(original_covariates[[j]])[!mask[, j]]
    if (length(unique(observed)) <= 2L) return(1)
    value <- stats::sd(observed)
    if (is.finite(value) && value > 1e-8) value else 1
  }, numeric(1))
  cell_scales <- unlist(lapply(seq_along(covariate_names), function(j) {
    rep(scale_by_column[j], sum(mask[, j]))
  }), use.names = FALSE)

  z <- sweep(draw_matrix, 2L, cell_scales, "/")
  z <- sweep(z, 2L, colMeans(z), "-")
  gram <- tcrossprod(z)
  norms <- diag(gram)
  d2 <- pmax(outer(norms, norms, "+") - 2 * gram, 0)
  distances <- sqrt(d2)
  total_distance <- rowSums(distances)
  local_id <- which.min(total_distance)
  draw_id <- keep[local_id]

  completed <- data
  selected <- fit$draws$X_mis[[draw_id]]
  completed[, covariate_names] <- selected[, covariate_names, drop = FALSE]
  list(
    data = completed,
    draw_id = draw_id,
    total_distance = total_distance[local_id],
    n_draws = nrow(draw_matrix),
    n_missing_cells = ncol(draw_matrix)
  )
}

.log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

.dinvgamma_log <- function(x, shape, rate) {
  shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x
}

# Chib's identity conditional on a fixed completed design matrix. The beta
# ordinate uses the full chain; the random-effects ordinate uses a reduced chain
# with beta fixed; the variance ordinate is an exact inverse-gamma density.
compute_chib_logml <- function(
  fit, 
  data, 
  theta_star,
  reduced_iter = 2000L,
  reduced_burnin = 500L,
  reduced_thin = 2L,
  seed = 202606L,
  details = FALSE
) {
  if (any(fit$draws$missing_mask)) {
    stop(
      "The posterior run used changing imputations. Select one common medoid, ",
      "refit this candidate model on that completed data set, and then call ",
      "compute_chib_logml().",
      call. = FALSE
    )
  }
  d <- posterior_draws_after_burnin(fit)
  design_data <- as.data.frame(data)
  design_data[[fit$id_var]] <- factor(design_data[[fit$id_var]])
  X <- model_matrix_from_formula(fit$fixed_formula, design_data)
  if (anyNA(X)) {
    stop("Chib evidence requires a fixed completed design matrix.", call. = FALSE)
  }

  y <- as.integer(design_data[[fit$y_var]])
  id <- as.integer(design_data[[fit$id_var]])
  J <- nlevels(design_data[[fit$id_var]])
  group_indices <- split(seq_along(id), id)
  beta_star <- as.numeric(theta_star$beta)
  a_star <- as.numeric(theta_star$a)
  sigma2_star <- as.numeric(theta_star$sigma2_a)
  stopifnot(length(beta_star) == ncol(X), length(a_star) == J)

  prior_precision <- fit$beta_prior$Sigma_inv
  prior_mean <- fit$beta_prior$mu
  V_beta <- solve(prior_precision + crossprod(X))
  log_beta_density <- vapply(seq_len(nrow(d$z)), function(s) {
    m_beta <- V_beta %*% (
      prior_precision %*% prior_mean +
        crossprod(X, d$z[s, ] - d$a[s, id])
    )
    mvtnorm::dmvnorm(beta_star, mean = as.numeric(m_beta),
                     sigma = V_beta, log = TRUE)
  }, numeric(1))
  log_ordinate_beta <- .log_mean_exp(log_beta_density)

  set.seed(seed)
  a <- a_star
  sigma2_a <- sigma2_star
  z_latent <- ifelse(y == 1L, 0.5, -0.5)
  Xbeta <- as.numeric(X %*% beta_star)
  log_a_density <- numeric(0)

  for (iter in seq_len(reduced_iter)) {
    eta <- Xbeta + a[id]
    z_latent[y == 1L] <- truncnorm::rtruncnorm(
      sum(y == 1L), a = 0, b = Inf, mean = eta[y == 1L], sd = 1
    )
    z_latent[y == 0L] <- truncnorm::rtruncnorm(
      sum(y == 0L), a = -Inf, b = 0, mean = eta[y == 0L], sd = 1
    )

    for (j in seq_len(J)) {
      idx <- group_indices[[j]]
      V_a <- 1 / (length(idx) + 1 / sigma2_a)
      m_a <- V_a * sum(z_latent[idx] - Xbeta[idx])
      a[j] <- stats::rnorm(1L, m_a, sqrt(V_a))
    }
    shape <- fit$c0 + J / 2
    rate <- fit$d0 + 0.5 * sum(a^2)
    sigma2_a <- 1 / stats::rgamma(1L, shape = shape, rate = rate)

    if (iter > reduced_burnin && (iter - reduced_burnin) %% reduced_thin == 0L) {
      value <- 0
      for (j in seq_len(J)) {
        idx <- group_indices[[j]]
        V_a <- 1 / (length(idx) + 1 / sigma2_a)
        m_a <- V_a * sum(z_latent[idx] - Xbeta[idx])
        value <- value + stats::dnorm(a_star[j], m_a, sqrt(V_a), log = TRUE)
      }
      log_a_density <- c(log_a_density, value)
    }
  }
  if (!length(log_a_density)) stop("Reduced Gibbs retained no draws.", call. = FALSE)
  log_ordinate_a <- .log_mean_exp(log_a_density)
  log_ordinate_sigma2 <- .dinvgamma_log(
    sigma2_star,
    shape = fit$c0 + J / 2,
    rate = fit$d0 + 0.5 * sum(a_star^2)
  )

  eta_star <- as.numeric(X %*% beta_star + a_star[id])
  prob_star <- pmin(pmax(stats::pnorm(eta_star), 1e-12), 1 - 1e-12)
  log_likelihood <- sum(stats::dbinom(y, 1L, prob_star, log = TRUE))
  log_prior <- mvtnorm::dmvnorm(
    beta_star, mean = prior_mean, sigma = fit$beta_prior$Sigma, log = TRUE
  ) + sum(stats::dnorm(a_star, 0, sqrt(sigma2_star), log = TRUE)) +
    .dinvgamma_log(sigma2_star, fit$c0, fit$d0)
  log_posterior_ordinate <- log_ordinate_beta + log_ordinate_a + log_ordinate_sigma2
  logml <- log_likelihood + log_prior - log_posterior_ordinate

  if (!isTRUE(details)) return(logml)
  list(
    log_marginal_likelihood = logml,
    log_likelihood_at_theta = log_likelihood,
    log_prior_at_theta = log_prior,
    log_posterior_ordinate = log_posterior_ordinate,
    log_ordinate_beta = log_ordinate_beta,
    log_ordinate_a_given_beta = log_ordinate_a,
    log_ordinate_sigma2_given_beta_a = log_ordinate_sigma2,
    reduced_draws = length(log_a_density),
    evidence_target = "p(y|X_medoid,M): conditional marginal likelihood"
  )
}

run_chib_fixed_point_sensitivity <- function(fit, data) {
  fps <- get_chib_fixed_points(fit)
  do.call(rbind, lapply(names(fps), function(nm) {
    logml <- tryCatch(compute_chib_logml(fit, data, fps[[nm]]), error = function(e) NA_real_)
    data.frame(fixed_point = nm, log_marginal_likelihood = logml)
  }))
}

run_prior_sensitivity <- function(data,
                                  fixed_formula,
                                  id_var = "id",
                                  y_var = "y",
                                  tau2_grid = c(0.5, 1, 2.5, 5, 10),
                                  mcmc_control,
                                  seed = 123) {
  res <- vector("list", length(tau2_grid))

  for (i in seq_along(tau2_grid)) {
    tau2 <- tau2_grid[i]
    tmp_design <- model_matrix_from_formula(fixed_formula, data)
    prior <- make_beta_prior(ncol(tmp_design), tau2_beta = tau2)

    timed <- with_runtime({
      fit_bayes_re_probit_cart(
        data = data,
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
        seed = seed + i,
        verbose_every = mcmc_control$verbose_every
      )
    })

    res[[i]] <- list(tau2_beta = tau2, fit = timed$value, runtime_seconds = timed$runtime_seconds)
  }

  res
}
