# -----------------------------------------------------------------------------
# 00_config.R
# Configuration for Probit models with missingness in covariates
# -----------------------------------------------------------------------------

required_packages <- c(
  "MASS", "truncnorm", "rpart", "lme4", "mice", "mitml", "mvtnorm", "coda"
)

optional_packages <- c("missForest", "loo")

load_required_packages <- function(pkgs = required_packages) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Install missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# run only once!
# Use relative paths by default. Override with options() or Sys.setenv().
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data")
output_dir  <- file.path(project_dir, "output")
results_dir <- file.path(output_dir, "results")
tables_dir  <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# Simulation defaults
sim_config <- list(
  n_rep = 100,
  N = 2000,
  J = 50,
  beta_true = c(
    "(Intercept)" = 0.20,
    "x2" = -1.50,
    "x3" = 0.80,
    "x4" = -0.30,
    "x5B" = -0.40,
    "x5C" = -0.30
  ),
  beta_full = c(
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
  ),
    
  sigma2_a_true = 1.00,
  missing_rates_mcar = c(
    x2  = 0.20,
    x3  = 0.30,
    x5B = 0.10,
    x5C = 0.15,
    x7  = 0.20,
    x10 = 0.30
  ),
  seed = 202605L
)


mcmc_control <- list(
  n_iter = 4000,
  burnin = 1000,
  thin = 1,
  verbose_every = 1000,
  cart_minbucket_range = 20:60,
  cart_maxdepth = 4,
  cart_cp = 0.001,
  cart_max_donors = 50
)

make_beta_prior <- function(p, tau2_beta = 1) {
  list(
    mu = rep(0, p),
    Sigma = diag(tau2_beta, p),
    Sigma_inv = diag(1 / tau2_beta, p),
    tau2_beta = tau2_beta
  )
}

sigma2_a_prior <- list(c0 = 2, d0 = 1)

# Candidate models for the simulated data. Adjust names/formulas if needed.
model_specs <- list(
  M1_intercept = list(fixed = y ~ 1, k_true = FALSE),
  M2_true      = list(fixed = y ~ x2 + x3 + x4 + x5B + x5C, k_true = TRUE),
  M3_reduced   = list(fixed = y ~ x2 + x5B + x5C, k_true = FALSE),
  M4_noise     = list(fixed = y ~ x4 + x5B + x5C + x7 + x8B + x8C, k_true = FALSE),
  M5_full      = list(fixed = y ~ x2 + x3 + x4 + x5B + x5C + x7 + x8B + x8C + x10, k_true = FALSE),
  M6_alt       = list(fixed = y ~ x2 + x5B + x5C + x7 + x10, k_true = FALSE)
)

# Prior variances grid
tau2_beta_grid <- c(0.5, 1, 2.5, 5, 10)
