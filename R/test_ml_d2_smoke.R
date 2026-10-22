# -----------------------------------------------------------------------------
# test_ml_d2_smoke.R
# Fast end-to-end check of the revised ML path. No Julia/Bayesian code is run.
# Run from the ZIP-upload project root with:
#   Rscript R/test_ml_d2_smoke.R
# or in RStudio with the project root as the working directory.
# -----------------------------------------------------------------------------

rm(list = ls())

source("R/00_config.R")
load_required_packages(c("MASS", "lme4", "mice", "mitml"))
source("R/01_simulation_dgp_missingness.R")
source("R/04_ml_benchmarks.R")
source("R/05_run_simulation.R")
source("R/08_ml_before_deletion.R")
source("R/09_ml_imputation_estimates.R")
source("R/10_mice_d2_lrt.R")

smoke_data <- make_simulation_data(
  n_rep = 1L,
  mechanisms = "MCAR"
)

smoke <- run_ml_benchmarks_one(
  sim_item = smoke_data[[1L]],
  model_specs = model_specs,
  run_missforest = FALSE,
  mice_m = 5L,
  mice_maxit = 2L,
  seed = 202605L
)

stopifnot(
  nrow(smoke$ml) > 0L,
  nrow(smoke$ml_estimates) > 0L,
  nrow(smoke$mice_d2_tests) == nrow(default_nested_model_pairs()),
  all(smoke$mice_d2_tests$pooling == "D2-likelihood"),
  any(is.finite(smoke$mice_d2_tests$p_value)),
  all(smoke$mice_d2_tests$n_imp >= 1L)
)

print(smoke$ml[, c("method", "model", "converged", "BIC")])
print(smoke$mice_d2_tests)
message("ML likelihood-based D2 smoke test passed.")
