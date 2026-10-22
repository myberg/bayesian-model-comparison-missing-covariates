###
# 73_SC6_postprocess_bayes_results.R
#
# Post-process Julia Bayesian CART-BB / common-medoid / Chib outputs.
# Produces coefficient tables, conditional log-BF matrices and compact
# model-selection summaries for the paper workflow.
###

message_ts <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

if (!exists("PROJECT_DIR", inherits = TRUE)) {
  PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SC6_BAYES_DIR <- Sys.getenv(
  "SC6_BAYES_DIR",
  unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "04_bayes_cartbb")
)
SC6_EVIDENCE_DIR <- Sys.getenv(
  "SC6_EVIDENCE_DIR",
  unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "05_conditional_evidence")
)

dir.create(SC6_EVIDENCE_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("dplyr", "tidyr", "readr", "tibble")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0L) {
  stop("Required package(s) missing: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

required_files <- file.path(
  SC6_BAYES_DIR,
  c(
    "bayes_model_results.csv",
    "bayes_estimates.csv",
    "bayes_marginal_likelihood.csv",
    "bayes_predictive_criteria.csv",
    "bayes_medoid_metadata.csv"
  )
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Julia Bayesian output file(s) missing:\n",
    paste(" -", missing_files, collapse = "\n"),
    "\nRun 71_SC6_run_bayes_medoid.jl first.",
    call. = FALSE
  )
}

model_results <- readr::read_csv(file.path(SC6_BAYES_DIR, "bayes_model_results.csv"), show_col_types = FALSE)
estimates <- readr::read_csv(file.path(SC6_BAYES_DIR, "bayes_estimates.csv"), show_col_types = FALSE)
marglik <- readr::read_csv(file.path(SC6_BAYES_DIR, "bayes_marginal_likelihood.csv"), show_col_types = FALSE)
predictive <- readr::read_csv(file.path(SC6_BAYES_DIR, "bayes_predictive_criteria.csv"), show_col_types = FALSE)
medoid_metadata <- readr::read_csv(file.path(SC6_BAYES_DIR, "bayes_medoid_metadata.csv"), show_col_types = FALSE)

candidate_models <- paste0("M", 1:5)
model_labels <- tibble::tribble(
  ~model, ~block,
  "M1", "Baseline",
  "M2", "Baseline + human capital",
  "M3", "Baseline + labour-market history",
  "M4", "Baseline + social background",
  "M5", "Baseline + all substantive blocks"
)

marglik_mean <- marglik %>%
  filter(.data$fixed_point == "mean") %>%
  arrange(match(.data$model, candidate_models))

if (nrow(marglik_mean) != length(candidate_models)) {
  stop("Expected one mean fixed-point Chib row for each of M1--M5.", call. = FALSE)
}

bayes_ranking <- marglik_mean %>%
  select(
    model, fixed_point, log_marginal_likelihood,
    log_likelihood_at_theta, log_prior_at_theta,
    log_posterior_ordinate, sigma2_a, nobs, n_groups,
    n_iter, burnin, thin
  ) %>%
  arrange(desc(.data$log_marginal_likelihood)) %>%
  mutate(
    rank_conditional_logml = row_number(),
    delta_logml_from_best = .data$log_marginal_likelihood - max(.data$log_marginal_likelihood),
    logBF_vs_best = .data$log_marginal_likelihood - max(.data$log_marginal_likelihood)
  ) %>%
  left_join(model_labels, by = "model") %>%
  select(model, block, everything())

readr::write_csv(bayes_ranking, file.path(SC6_EVIDENCE_DIR, "bayes_conditional_logml_ranking.csv"))

logml_vec <- marglik_mean$log_marginal_likelihood
names(logml_vec) <- marglik_mean$model
logml_vec <- logml_vec[candidate_models]

logbf_mat <- outer(logml_vec, logml_vec, FUN = "-")
rownames(logbf_mat) <- candidate_models
colnames(logbf_mat) <- candidate_models

logbf_matrix <- as.data.frame(logbf_mat) %>%
  tibble::rownames_to_column("model_i")

logbf_pairwise <- expand.grid(
  model_i = candidate_models,
  model_j = candidate_models,
  stringsAsFactors = FALSE
) %>%
  mutate(logBF_i_vs_j = logbf_mat[cbind(.data$model_i, .data$model_j)])

readr::write_csv(logbf_matrix, file.path(SC6_EVIDENCE_DIR, "bayes_conditional_logbf_matrix.csv"))
readr::write_csv(logbf_pairwise, file.path(SC6_EVIDENCE_DIR, "bayes_conditional_logbf_pairwise.csv"))

coef_wide <- estimates %>%
  filter(.data$component %in% c("fixed", "random_effect_variance")) %>%
  mutate(
    term = dplyr::recode(
      .data$term,
      "(Intercept)" = "Intercept",
      "sigma2_a" = "Random-intercept variance",
      "unemp_lt1" = "Unemployment duration: <1 year",
      "unemp_1to2" = "Unemployment duration: 1 to <2 years",
      "unemp_2to5" = "Unemployment duration: 2 to <5 years",
      "unemp_ge5" = "Unemployment duration: >=5 years",
      .default = .data$term
    ),
    estimate_ci = sprintf("%.3f [%.3f, %.3f]", .data$posterior_mean, .data$ci_lower, .data$ci_upper)
  ) %>%
  select(term, model, estimate_ci) %>%
  tidyr::pivot_wider(names_from = .data$model, values_from = .data$estimate_ci) %>%
  arrange(.data$term)

readr::write_csv(coef_wide, file.path(SC6_EVIDENCE_DIR, "paper_bayes_coefficients_wide.csv"))

predictive_summary <- predictive %>%
  filter(.data$model %in% candidate_models) %>%
  select(model, WAIC, LOO, lppd, p_waic, elpd_waic, elpd_loo, p_loo, n_draws) %>%
  arrange(.data$WAIC) %>%
  mutate(rank_WAIC = row_number()) %>%
  left_join(model_labels, by = "model") %>%
  select(model, block, everything())

readr::write_csv(predictive_summary, file.path(SC6_EVIDENCE_DIR, "paper_bayes_predictive_criteria.csv"))

# Fixed-point sensitivity, if Julia was run with several fixed points.
fixed_point_sensitivity <- marglik %>%
  select(model, fixed_point, log_marginal_likelihood, sigma2_a, reduced_draws) %>%
  arrange(.data$model, .data$fixed_point)
readr::write_csv(fixed_point_sensitivity, file.path(SC6_EVIDENCE_DIR, "bayes_fixed_point_sensitivity.csv"))

# Lightweight LaTeX helper
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
  rows <- apply(df_chr, 1, function(z) paste(z, collapse = " & "))
  lines <- c(
    "\\begin{table}[tbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\small",
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0(paste(escape_latex(names(df_chr)), collapse = " & "), " \\\\"),
    "\\midrule",
    paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, path)
}

bayes_rank_tex <- bayes_ranking %>%
  transmute(
    Model = .data$model,
    Block = .data$block,
    `log p(y|X_med,M)` = sprintf("%.2f", .data$log_marginal_likelihood),
    `Rank` = .data$rank_conditional_logml,
    `log BF vs best` = sprintf("%.2f", .data$logBF_vs_best)
  )

write_simple_latex_table(
  bayes_rank_tex,
  file.path(SC6_EVIDENCE_DIR, "table_bayes_conditional_logml_ranking.tex"),
  caption = "NEPS SC6 Bayesian workflow: conditional Chib evidence on the common CART-BB medoid.",
  label = "tab:sc6_bayes_conditional_logml"
)

message_ts("Bayesian post-processing written to: ", SC6_EVIDENCE_DIR)
message_ts("Primary files:")
message_ts(" - bayes_conditional_logml_ranking.csv")
message_ts(" - bayes_conditional_logbf_matrix.csv")
message_ts(" - bayes_conditional_logbf_pairwise.csv")
message_ts(" - paper_bayes_coefficients_wide.csv")
message_ts(" - paper_bayes_predictive_criteria.csv")
