# -----------------------------------------------------------------------------
# 16_make_traceplots.R
# Traceplot generation for selected Bayesian workflows.
#
# Expected input:
#   outputs/<run_id>/06_mcmc_diagnostics/draws_*.csv
#
# Output:
#   outputs/<run_id>/08_figures/traceplots_*.pdf
#   outputs/<run_id>/08_figures_data/plotdata_traceplots.csv
# -----------------------------------------------------------------------------

make_traceplots <- function(root,
                            terms = c("(Intercept)", "x2", "x3", "x4", "x5B", "x5C", "sigma2_a", "loglik"),
                            workflows = NULL,
                            mechanisms = c("MCAR", "MAR"),
                            models = c("M2_true"),
                            width = 10,
                            height = 7) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for traceplots.")
  }

  diag_dir <- file.path(root, "06_mcmc_diagnostics")
  fig_dir <- file.path(root, "08_figures")
  fig_data_dir <- file.path(root, "08_figures_data")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_data_dir, recursive = TRUE, showWarnings = FALSE)

  if (!exists("read_mcmc_draws", mode = "function")) {
    source("R/15_make_mcmc_diagnostics.R")
  }

  draws <- read_mcmc_draws(root)
  if (is.null(draws)) return(invisible(NULL))

  plotdata <- draws |>
    dplyr::filter(.data$term %in% terms,
                  .data$mechanism %in% mechanisms,
                  .data$model %in% models)

  if (!is.null(workflows)) {
    plotdata <- plotdata |>
      dplyr::filter(.data$workflow %in% workflows)
  }

  if (nrow(plotdata) == 0L) {
    warning("No traceplot data after filtering.")
    return(invisible(NULL))
  }

  readr::write_csv(plotdata, file.path(fig_data_dir, "plotdata_traceplots.csv"))

  groups <- plotdata |>
    dplyr::distinct(.data$workflow, .data$mechanism, .data$model, .data$rep)

  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    dat_g <- plotdata |>
      dplyr::filter(.data$workflow == g$workflow,
                    .data$mechanism == g$mechanism,
                    .data$model == g$model,
                    .data$rep == g$rep)

    p <- ggplot2::ggplot(dat_g, ggplot2::aes(x = .data$iter, y = .data$value)) +
      ggplot2::geom_line(linewidth = 0.25) +
      ggplot2::facet_wrap(~ term, scales = "free_y") +
      ggplot2::labs(
        x = "MCMC iteration (retained draws)",
        y = "Draw value",
        title = paste("Traceplots:", g$workflow, g$mechanism, g$model, paste0("rep=", g$rep))
      ) +
      ggplot2::theme_bw()

    file_stub <- paste(
      "traceplots",
      gsub("[^A-Za-z0-9]+", "_", tolower(g$workflow)),
      tolower(g$mechanism),
      g$model,
      sprintf("rep%03d", as.integer(g$rep)),
      sep = "_"
    )

    ggplot2::ggsave(
      filename = file.path(fig_dir, paste0(file_stub, ".pdf")),
      plot = p,
      width = width,
      height = height
    )
  }

  message("Traceplots written to: ", fig_dir)
  invisible(plotdata)
}
