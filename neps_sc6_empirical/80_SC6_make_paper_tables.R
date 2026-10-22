###
# 80_SC6_make_paper_tables.R
#
# Combine ML and Bayesian model-selection results into paper-ready overview
# tables for the NEPS SC6 Wave-16 empirical illustration.
###

message_ts <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

if (!exists("PROJECT_DIR", inherits = TRUE)) {
  PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SC6_ML_DIR <- Sys.getenv(
  "SC6_ML_DIR",
  unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "03_ml_benchmark")
)
SC6_EVIDENCE_DIR <- Sys.getenv(
  "SC6_EVIDENCE_DIR",
  unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "05_conditional_evidence")
)
SC6_TABLE_DIR <- Sys.getenv(
  "SC6_TABLE_DIR",
  unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "06_paper_tables")
)

dir.create(SC6_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

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

need <- c(
  file.path(SC6_ML_DIR, "ml_fit_statistics_summary.csv"),
  file.path(SC6_ML_DIR, "ml_bic_model_ranking.csv"),
  file.path(SC6_ML_DIR, "ml_bic_logbf_matrix.csv"),
  file.path(SC6_ML_DIR, "ml_lrt_nested_pooled.csv"),
  file.path(SC6_EVIDENCE_DIR, "bayes_conditional_logml_ranking.csv"),
  file.path(SC6_EVIDENCE_DIR, "bayes_conditional_logbf_matrix.csv"),
  file.path(SC6_EVIDENCE_DIR, "paper_bayes_predictive_criteria.csv")
)
missing <- need[!file.exists(need)]
if (length(missing) > 0L) {
  stop("Required model-comparison file(s) missing:\n", paste(" -", missing, collapse = "\n"), call. = FALSE)
}

ml_fit <- readr::read_csv(file.path(SC6_ML_DIR, "ml_fit_statistics_summary.csv"), show_col_types = FALSE)
ml_rank <- readr::read_csv(file.path(SC6_ML_DIR, "ml_bic_model_ranking.csv"), show_col_types = FALSE)
ml_bf <- readr::read_csv(file.path(SC6_ML_DIR, "ml_bic_logbf_matrix.csv"), show_col_types = FALSE)
ml_lrt <- readr::read_csv(file.path(SC6_ML_DIR, "ml_lrt_nested_pooled.csv"), show_col_types = FALSE)

bayes_rank <- readr::read_csv(file.path(SC6_EVIDENCE_DIR, "bayes_conditional_logml_ranking.csv"), show_col_types = FALSE)
bayes_bf <- readr::read_csv(file.path(SC6_EVIDENCE_DIR, "bayes_conditional_logbf_matrix.csv"), show_col_types = FALSE)
bayes_pred <- readr::read_csv(file.path(SC6_EVIDENCE_DIR, "paper_bayes_predictive_criteria.csv"), show_col_types = FALSE)

candidate_models <- paste0("M", 1:5)
model_labels <- tibble::tribble(
  ~model, ~block,
  "M1", "Baseline",
  "M2", "Baseline + human capital",
  "M3", "Baseline + labour-market history",
  "M4", "Baseline + social background",
  "M5", "Baseline + all substantive blocks"
)

ml_overview <- ml_fit %>%
  select(model, nobs, logLik_median, AIC_median, BIC_median, re_var_median) %>%
  left_join(ml_rank %>% select(model, rank_BIC), by = "model") %>%
  mutate(
    ML_delta_BIC = .data$BIC_median - min(.data$BIC_median),
    ML_logBF_vs_best = -0.5 * .data$ML_delta_BIC
  )

bayes_overview <- bayes_rank %>%
  select(model, Bayes_rank_logml = rank_conditional_logml,
         Bayes_logml = log_marginal_likelihood,
         Bayes_delta_logml = delta_logml_from_best,
         Bayes_logBF_vs_best = logBF_vs_best,
         Bayes_sigma2_a = sigma2_a) %>%
  left_join(
    bayes_pred %>% select(model, WAIC, rank_WAIC),
    by = "model"
  )

overview <- model_labels %>%
  left_join(ml_overview, by = "model") %>%
  left_join(bayes_overview, by = "model") %>%
  arrange(match(.data$model, candidate_models))

readr::write_csv(overview, file.path(SC6_TABLE_DIR, "paper_model_selection_overview.csv"))

# Pairwise logBF matrices: reshape to one row per model pair.
ml_pair <- ml_bf %>%
  rename(model_i = .data$model_i) %>%
  pivot_longer(all_of(candidate_models), names_to = "model_j", values_to = "ML_BIC_logBF_i_vs_j")

bayes_pair <- bayes_bf %>%
  rename(model_i = .data$model_i) %>%
  pivot_longer(all_of(candidate_models), names_to = "model_j", values_to = "Bayes_conditional_logBF_i_vs_j")

pairwise <- ml_pair %>%
  left_join(bayes_pair, by = c("model_i", "model_j")) %>%
  filter(.data$model_i != .data$model_j) %>%
  arrange(match(.data$model_i, candidate_models), match(.data$model_j, candidate_models))

readr::write_csv(pairwise, file.path(SC6_TABLE_DIR, "paper_pairwise_logbf_ml_bayes.csv"))

# Nested comparison table for main/appendix reporting.
lrt_compact <- ml_lrt %>%
  transmute(
    small, large, comparison,
    ML_LR_median = .data$LR_median,
    ML_df = .data$df_diff_median,
    ML_p_median = .data$p_median,
    ML_D2_p = .data$D2_p_value,
    ML_D2_status = .data$D2_status
  )
readr::write_csv(lrt_compact, file.path(SC6_TABLE_DIR, "paper_ml_nested_lrt_compact.csv"))

# LaTeX helper ------------------------------------------------------------

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("<", "$<$", x, fixed = TRUE)
  x <- gsub(">", "$>$", x, fixed = TRUE)
  x
}

write_simple_latex_table <- function(df, path, caption, label, align = NULL) {
  df_chr <- df %>% mutate(across(everything(), as.character))
  df_chr[] <- lapply(df_chr, escape_latex)
  if (is.null(align)) align <- paste0("l", paste(rep("r", ncol(df_chr) - 1L), collapse = ""))
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

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

overview_tex <- overview %>%
  transmute(
    Model = .data$model,
    Block = .data$block,
    `ML BIC` = fmt(.data$BIC_median, 2),
    `ML rank` = .data$rank_BIC,
    `ML log BF vs best` = fmt(.data$ML_logBF_vs_best, 2),
    `Bayes log evidence` = fmt(.data$Bayes_logml, 2),
    `Bayes rank` = .data$Bayes_rank_logml,
    `Bayes log BF vs best` = fmt(.data$Bayes_logBF_vs_best, 2),
    WAIC = fmt(.data$WAIC, 2)
  )

write_simple_latex_table(
  overview_tex,
  file.path(SC6_TABLE_DIR, "table_model_selection_overview.tex"),
  caption = "NEPS SC6 model-selection overview: ML benchmark and Bayesian conditional evidence.",
  label = "tab:sc6_model_selection_overview",
  align = "llrrrrrrr"
)

# A concise text summary with the selected models.
best_ml <- overview %>% arrange(.data$rank_BIC) %>% slice(1)
best_bayes <- overview %>% arrange(.data$Bayes_rank_logml) %>% slice(1)

summary_lines <- c(
  paste0("Best ML model by median BIC: ", best_ml$model, " (", best_ml$block, ")."),
  paste0("Best Bayesian model by conditional Chib evidence: ", best_bayes$model, " (", best_bayes$block, ")."),
  "Positive log BF values in pairwise tables favour model_i over model_j.",
  "Bayesian log evidence is conditional on the common CART-BB medoid completion."
)
writeLines(summary_lines, con = file.path(SC6_TABLE_DIR, "paper_model_selection_summary.txt"))

message_ts("Paper model-selection tables written to: ", SC6_TABLE_DIR)
message_ts("Primary files:")
message_ts(" - paper_model_selection_overview.csv")
message_ts(" - paper_pairwise_logbf_ml_bayes.csv")
message_ts(" - table_model_selection_overview.tex")
message_ts(" - paper_model_selection_summary.txt")
