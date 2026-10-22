###
# 91_SC6_bayes_sensitivity_driver_v4.R
#
# Bayesian sensitivity analyses for the empirical NEPS SC6 illustration.
# It reruns the Julia common-medoid workflow for a small, defensible grid of
# prior and medoid-reference checks and summarizes the conditional evidence
# ranking. Keep this separate from the main paper run.
###

message_ts <- function(...) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)

if (!exists("PROJECT_DIR", inherits = TRUE)) {
  PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SC6_APPENDIX_DIR <- Sys.getenv("SC6_APPENDIX_DIR", unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "07_appendix"))
review_dir <- file.path(SC6_APPENDIX_DIR, "reviewer_quality")
dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("dplyr", "readr", "tibble", "tidyr")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0L) {
  stop("Required package(s) missing: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

julia_bin <- Sys.getenv("JULIA_BIN", unset = Sys.which("julia"))
if (julia_bin == "") stop("Julia was not found. Set JULIA_BIN or add Julia to PATH.", call. = FALSE)

runner_rel <- file.path("neps_sc6_empirical", "71_SC6_run_bayes_medoid.jl")
runner_abs <- file.path(PROJECT_DIR, runner_rel)
if (!file.exists(runner_abs)) stop("Julia runner not found: ", runner_abs, call. = FALSE)

input_files <- list.files(PROJECT_DIR, pattern = "sc6_wave16_bayes_input\\.csv$", recursive = TRUE, full.names = FALSE)
if (length(input_files) == 0L) {
  stop("No sc6_wave16_bayes_input.csv found. Run 70_SC6_bayes_cartbb_models.R first.", call. = FALSE)
}
input_info <- file.info(file.path(PROJECT_DIR, input_files))
input_rel <- input_files[which.max(input_info$mtime)]
message_ts("Using input: ", input_rel)

# Main sensitivity grid. This is intentionally small: it checks whether the
# evidence ranking is stable to one tighter prior, one weaker prior and a
# different medoid-reference model. Expand only if rankings change.
sensitivity_grid <- tibble::tribble(
  ~scenario, ~tau2_beta, ~medoid_reference_model, ~seed,
  "main_ref_M5_tau1", 1.0, "M5", 20260628L,
  "prior_tighter_tau0p5", 0.5, "M5", 20260628L,
  "prior_weaker_tau2p5", 2.5, "M5", 20260628L,
  "medoid_ref_M3_tau1", 1.0, "M3", 20260628L
)

# Defaults suitable for a publishable appendix sensitivity check. For a quick
# technical test, temporarily reduce these values but do not report quick-trial
# results in the manuscript.
controls <- list(
  medoid_n_iter = as.integer(Sys.getenv("SC6_SENS_MEDOID_N_ITER", unset = "12000")),
  medoid_burnin = as.integer(Sys.getenv("SC6_SENS_MEDOID_BURNIN", unset = "3000")),
  medoid_thin = as.integer(Sys.getenv("SC6_SENS_MEDOID_THIN", unset = "10")),
  n_iter = as.integer(Sys.getenv("SC6_SENS_N_ITER", unset = "20000")),
  burnin = as.integer(Sys.getenv("SC6_SENS_BURNIN", unset = "5000")),
  thin = as.integer(Sys.getenv("SC6_SENS_THIN", unset = "5")),
  chib_reduced_iter = as.integer(Sys.getenv("SC6_SENS_CHIB_ITER", unset = "8000")),
  chib_reduced_burnin = as.integer(Sys.getenv("SC6_SENS_CHIB_BURNIN", unset = "2000")),
  chib_reduced_thin = as.integer(Sys.getenv("SC6_SENS_CHIB_THIN", unset = "2")),
  cart_every = 10L,
  cart_max_depth = 3L,
  cart_min_samples_leaf = 20L,
  cart_max_donors = 50L,
  fixed_points = "mean",
  verbose_every = 1000L
)

force_rerun <- tolower(Sys.getenv("SC6_SENS_FORCE_RERUN", unset = "false")) %in% c("true", "1", "yes")

run_one <- function(row) {
  out_rel <- file.path("results", "sc6_empirical", "bayes_sensitivity", row$scenario)
  out_abs <- file.path(PROJECT_DIR, out_rel)
  dir.create(out_abs, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    "--project=.",
    runner_rel,
    "--data", input_rel,
    "--outdir", out_rel,
    "--medoid_reference_model", row$medoid_reference_model,
    "--tau2_beta", as.character(row$tau2_beta),
    "--seed", as.character(row$seed),
    "--medoid_n_iter", as.character(controls$medoid_n_iter),
    "--medoid_burnin", as.character(controls$medoid_burnin),
    "--medoid_thin", as.character(controls$medoid_thin),
    "--n_iter", as.character(controls$n_iter),
    "--burnin", as.character(controls$burnin),
    "--thin", as.character(controls$thin),
    "--chib_reduced_iter", as.character(controls$chib_reduced_iter),
    "--chib_reduced_burnin", as.character(controls$chib_reduced_burnin),
    "--chib_reduced_thin", as.character(controls$chib_reduced_thin),
    "--cart_every", as.character(controls$cart_every),
    "--cart_max_depth", as.character(controls$cart_max_depth),
    "--cart_min_samples_leaf", as.character(controls$cart_min_samples_leaf),
    "--cart_max_donors", as.character(controls$cart_max_donors),
    "--fixed_points", controls$fixed_points,
    "--verbose_every", as.character(controls$verbose_every)
  )

  ranking_file <- file.path(out_abs, "bayes_conditional_logml_ranking.csv")
  medoid_file <- file.path(out_abs, "bayes_medoid_metadata.csv")

  if (!file.exists(ranking_file) || isTRUE(force_rerun)) {
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(PROJECT_DIR)
    message_ts("Running sensitivity scenario: ", row$scenario)
    status <- system2(command = julia_bin, args = args, stdout = "", stderr = "")
    if (!identical(status, 0L)) warning("Julia exit status for ", row$scenario, ": ", status, call. = FALSE)
  } else {
    message_ts("Using existing sensitivity output for scenario: ", row$scenario)
  }
  if (!file.exists(ranking_file)) {
    return(tibble(scenario = row$scenario, tau2_beta = row$tau2_beta, medoid_reference_model = row$medoid_reference_model,
                  seed = row$seed, status = "failed_no_ranking"))
  }

  ranking_df <- readr::read_csv(ranking_file, show_col_types = FALSE)
  medoid <- if (file.exists(medoid_file)) readr::read_csv(medoid_file, show_col_types = FALSE) else tibble()

  required_cols <- c("model", "log_marginal_likelihood")
  missing_cols <- setdiff(required_cols, names(ranking_df))
  if (length(missing_cols) > 0L) {
    return(tibble(
      scenario = row$scenario, tau2_beta = row$tau2_beta,
      medoid_reference_model = row$medoid_reference_model, seed = row$seed,
      status = paste0("failed_missing_columns: ", paste(missing_cols, collapse = ","))
    ))
  }

  if (!"rank_conditional_logml" %in% names(ranking_df)) {
    ranking_df <- ranking_df %>% arrange(desc(.data$log_marginal_likelihood)) %>%
      mutate(rank_conditional_logml = row_number())
  }

  ranking_df <- ranking_df %>% arrange(.data$rank_conditional_logml)
  best <- ranking_df %>% slice(1)
  logml_M3 <- ranking_df$log_marginal_likelihood[ranking_df$model == "M3"]
  logml_M5 <- ranking_df$log_marginal_likelihood[ranking_df$model == "M5"]
  m3_rank <- ranking_df$rank_conditional_logml[ranking_df$model == "M3"]
  m5_rank <- ranking_df$rank_conditional_logml[ranking_df$model == "M5"]
  ranking_string <- paste(ranking_df$model, collapse = ", ")

  tibble(
    scenario = row$scenario,
    tau2_beta = row$tau2_beta,
    medoid_reference_model = row$medoid_reference_model,
    seed = row$seed,
    status = "ok",
    medoid_draw_id = if (nrow(medoid) > 0 && "medoid_draw_id" %in% names(medoid)) medoid$medoid_draw_id[1] else NA_integer_,
    total_scaled_distance = if (nrow(medoid) > 0 && "total_scaled_distance" %in% names(medoid)) medoid$total_scaled_distance[1] else NA_real_,
    best_model = best$model[1],
    ranking = ranking_string,
    logml_best = best$log_marginal_likelihood[1],
    logBF_M3_vs_M5 = if (length(logml_M3) == 1L && length(logml_M5) == 1L) logml_M3 - logml_M5 else NA_real_,
    M3_rank = if (length(m3_rank) == 1L) m3_rank else NA_integer_,
    M5_rank = if (length(m5_rank) == 1L) m5_rank else NA_integer_
  )
}

summary <- bind_rows(lapply(seq_len(nrow(sensitivity_grid)), function(i) run_one(sensitivity_grid[i, ])))
readr::write_csv(summary, file.path(review_dir, "bayes_prior_medoid_sensitivity_summary.csv"))

fmt <- function(x, digits = 2) ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
escape_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x
}
tex <- summary %>%
  transmute(
    Scenario = scenario,
    `$\\tau_\\beta^2$` = fmt(tau2_beta, 1),
    Reference = medoid_reference_model,
    Medoid = medoid_draw_id,
    Best = best_model,
    Ranking = ranking,
    `$\\log BF_{3,5}$` = fmt(logBF_M3_vs_M5, 2)
  )
tex_chr <- tex %>% mutate(across(everything(), escape_latex))
rows <- if (nrow(tex_chr) == 0L) character(0) else apply(tex_chr, 1, function(z) paste(z, collapse = " & "))
lines <- c(
  "\\begin{table}[tbp]",
  "\\centering",
  "\\begin{threeparttable}",
  "\\caption{Prior and medoid-reference sensitivity of Bayesian conditional evidence.}",
  "\\label{tab:sc6_bayes_prior_medoid_sensitivity}",
  "\\small",
  "\\begin{tabular}{@{}llrlllr@{}}",
  "\\toprule",
  paste0(paste(names(tex_chr), collapse = " & "), " \\\\"),
  "\\midrule",
  paste0(rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  "\\item Notes: Positive values of $\\log BF_{3,5}$ favour the labour-market-history model over the encompassing model. All entries are conditional on the scenario-specific CART--BB medoid completion.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(lines, file.path(review_dir, "bayes_prior_medoid_sensitivity_summary.tex"))

message_ts("Bayesian sensitivity summary written to: ", normalizePath(review_dir, winslash = "/", mustWork = FALSE))
