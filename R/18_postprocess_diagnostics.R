# -----------------------------------------------------------------------------
# 18_postprocess_diagnostics.R
# One-call post-processing wrapper for reviewer-relevant diagnostics.
# -----------------------------------------------------------------------------

postprocess_reviewer_diagnostics <- function(root,
                                             make_plots = TRUE,
                                             trace_terms = c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a", "loglik")) {
  source("R/15_make_mcmc_diagnostics.R")
  source("R/16_make_traceplots.R")
  source("R/17_sensitivity_summaries.R")

  make_fixed_point_sensitivity(root)
  make_prior_sensitivity_tau2(root)  # emits warning if no prior_evidence_by_replication.csv exists yet
  make_mcmc_diagnostics(root)        # emits warning if no draws_*.csv exists yet

  if (isTRUE(make_plots)) {
    make_traceplots(root, terms = trace_terms)
  }

  invisible(TRUE)
}
