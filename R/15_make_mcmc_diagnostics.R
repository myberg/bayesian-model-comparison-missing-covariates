# -----------------------------------------------------------------------------
# 15_make_mcmc_diagnostics.R
# MCMC diagnostics for the SMAP/SMA workflow.
#
# Expected input:
#   outputs/<run_id>/06_mcmc_diagnostics/draws_*.csv
#
# Expected draw-file columns:
#   iter, chain, rep, mechanism, workflow, model, term, value
#
# Output:
#   ess_summary.csv
#   autocorr_summary.csv
#   mcse_summary.csv
#   diagnostics_summary.csv
# -----------------------------------------------------------------------------

safe_autocorr <- function(x, lag = 1L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= lag + 2L) return(NA_real_)
  stats::cor(x[seq_len(n - lag)], x[(lag + 1L):n], use = "complete.obs")
}

# A simple single-chain initial-positive-sequence ESS approximation.
# This is intentionally dependency-light and suitable for smoke-test diagnostics.
# If multiple chains are available, compute this per chain or adapt to posterior/coda.
approx_ess_one_chain <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 10L) return(NA_real_)
  if (stats::sd(x) == 0 || !is.finite(stats::sd(x))) return(NA_real_)

  if (is.null(max_lag)) {
    max_lag <- min(1000L, floor(n / 2))
  }

  ac <- vapply(seq_len(max_lag), function(k) safe_autocorr(x, lag = k), numeric(1))
  ac <- ac[is.finite(ac)]
  if (length(ac) == 0L) return(n)

  # Initial positive sequence: stop after the first non-positive autocorrelation.
  pos <- ac[ac > 0]
  if (length(pos) == 0L) return(n)

  # More conservative: use only consecutive positive values from the start.
  first_nonpos <- which(ac <= 0)[1]
  if (!is.na(first_nonpos) && first_nonpos > 1L) {
    pos <- ac[seq_len(first_nonpos - 1L)]
  } else if (!is.na(first_nonpos) && first_nonpos == 1L) {
    pos <- numeric(0)
  }

  tau <- 1 + 2 * sum(pos, na.rm = TRUE)
  ess <- n / tau
  max(1, min(n, ess))
}

compute_diag_one_series <- function(df, max_lag = 50L) {
  x <- df$value
  ess <- approx_ess_one_chain(x)
  sd_x <- stats::sd(x, na.rm = TRUE)
  mcse <- if (is.finite(ess) && ess > 0) sd_x / sqrt(ess) else NA_real_

  tibble::tibble(
    n_draws = sum(is.finite(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd_x,
    q025 = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q50 = stats::quantile(x, 0.500, na.rm = TRUE, names = FALSE),
    q975 = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE),
    ess_approx = ess,
    mcse_mean = mcse,
    autocorr_lag1 = safe_autocorr(x, 1L),
    autocorr_lag5 = safe_autocorr(x, 5L),
    autocorr_lag10 = safe_autocorr(x, 10L)
  )
}

read_mcmc_draws <- function(root) {
  diag_dir <- file.path(root, "06_mcmc_diagnostics")
  files <- list.files(diag_dir, pattern = "^draws_.*\\.csv$", full.names = TRUE)

  if (length(files) == 0L) {
    warning("No MCMC draw files found in ", diag_dir)
    return(NULL)
  }

  draws <- lapply(files, readr::read_csv, show_col_types = FALSE) |>
    dplyr::bind_rows()

  required <- c("iter", "chain", "rep", "mechanism", "workflow", "model", "term", "value")
  missing_cols <- setdiff(required, names(draws))
  if (length(missing_cols) > 0L) {
    stop("Draw files are missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  draws |>
    dplyr::mutate(
      mechanism = toupper(.data$mechanism),
      workflow = as.character(.data$workflow),
      model = as.character(.data$model),
      term = as.character(.data$term),
      value = as.numeric(.data$value)
    )
}

make_mcmc_diagnostics <- function(root) {
  diag_dir <- file.path(root, "06_mcmc_diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  draws <- read_mcmc_draws(root)
  if (is.null(draws)) return(invisible(NULL))

  diagnostics <- draws |>
    dplyr::group_by(.data$workflow, .data$mechanism, .data$model, .data$rep, .data$chain, .data$term) |>
    dplyr::group_modify(~ compute_diag_one_series(.x)) |>
    dplyr::ungroup()

  ess_summary <- diagnostics |>
    dplyr::select(.data$workflow, .data$mechanism, .data$model, .data$rep, .data$chain,
                  .data$term, .data$n_draws, .data$ess_approx)

  autocorr_summary <- diagnostics |>
    dplyr::select(.data$workflow, .data$mechanism, .data$model, .data$rep, .data$chain,
                  .data$term, .data$autocorr_lag1, .data$autocorr_lag5, .data$autocorr_lag10)

  mcse_summary <- diagnostics |>
    dplyr::select(.data$workflow, .data$mechanism, .data$model, .data$rep, .data$chain,
                  .data$term, .data$mcse_mean)

  readr::write_csv(ess_summary, file.path(diag_dir, "ess_summary.csv"))
  readr::write_csv(autocorr_summary, file.path(diag_dir, "autocorr_summary.csv"))
  readr::write_csv(mcse_summary, file.path(diag_dir, "mcse_summary.csv"))
  readr::write_csv(diagnostics, file.path(diag_dir, "diagnostics_summary.csv"))

  message("MCMC diagnostics written to: ", diag_dir)
  invisible(list(
    draws = draws,
    diagnostics = diagnostics,
    ess_summary = ess_summary,
    autocorr_summary = autocorr_summary,
    mcse_summary = mcse_summary
  ))
}
