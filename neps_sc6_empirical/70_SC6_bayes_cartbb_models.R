###
# 70_SC6_bayes_cartbb_models.R
#
# Export the prepared NEPS SC6 Wave-16 analysis data to a numeric Julia input
# file and optionally launch the Julia CART-BB common-medoid Bayesian workflow.
#
# This script is intentionally light on R-side modelling.  The Bayesian sampler,
# CART-BB donor imputation, common-medoid selection and conditional Chib evidence
# are handled in Julia by 71_SC6_run_bayes_medoid.jl and src/BayesREProbit.jl.
#
# To run Julia automatically from this R stage, set:
#   Sys.setenv(SC6_RUN_JULIA = "true")
# Otherwise, the script writes ready-to-run .bat/.sh command files.
###

message_ts <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

if (!exists("PROJECT_DIR", inherits = TRUE)) {
  PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
if (!exists("SC6_CONTROL", inherits = TRUE)) {
  SC6_CONTROL <- list(
    outcome = "curr_empl",
    group = "state_id",
    bayes_n_iter = 20000L,
    bayes_burnin = 4000L,
    bayes_thin = 2L,
    cart_every = 10L,
    cart_max_depth = 5L,
    cart_min_samples_leaf = 10L,
    cart_max_donors = 50L,
    medoid_reference_model = "M5",
    medoid_n_iter = 20000L,
    medoid_burnin = 5000L,
    medoid_thin = 10L,
    chib_reduced_iter = 4000L,
    chib_reduced_burnin = 1000L,
    chib_reduced_thin = 2L,
    tau2_beta = 1.0,
    c0 = 2.0,
    d0 = 1.0,
    seed = 20260628L
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

SC6_OUTPUT_DIR <- Sys.getenv(
  "SC6_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "data", "neps_sc6", "processed")
)
SC6_BAYES_DIR <- Sys.getenv(
  "SC6_BAYES_DIR",
  unset = file.path(PROJECT_DIR, "results", "sc6_empirical", "manual_run", "04_bayes_cartbb")
)
SC6_SCRIPT_DIR <- Sys.getenv("SC6_SCRIPT_DIR", unset = PROJECT_DIR)

dir.create(SC6_BAYES_DIR, recursive = TRUE, showWarnings = FALSE)
input_dir <- file.path(SC6_BAYES_DIR, "00_input")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("dplyr", "readr", "tibble")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0L) {
  stop("Required package(s) missing: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

analysis_file <- file.path(SC6_OUTPUT_DIR, "SC6_wave16_employment_analysis.rds")
model_variables_file <- file.path(SC6_OUTPUT_DIR, "SC6_wave16_model_variables.rds")

if (!file.exists(analysis_file)) {
  stop("Prepared Wave-16 analysis file not found: ", analysis_file, call. = FALSE)
}
if (!file.exists(model_variables_file)) {
  stop("Model-variable file not found: ", model_variables_file, call. = FALSE)
}

dat <- readRDS(analysis_file)
model_variables <- readRDS(model_variables_file)

required_cols <- unique(c(
  "ID_t", "curr_empl", "state_id", "unemployment_duration_cat",
  unlist(model_variables[paste0("M", 1:5)], use.names = FALSE)
))
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0L) {
  stop("Prepared data are missing required Bayesian columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

unemp_levels <- c("no previous spell", "<1 year", "1 to <2 years", "2 to <5 years", ">=5 years")
if (!all(na.omit(unique(as.character(dat$unemployment_duration_cat))) %in% unemp_levels)) {
  observed <- sort(unique(as.character(dat$unemployment_duration_cat)))
  stop(
    "Unexpected unemployment_duration_cat level(s): ",
    paste(setdiff(observed, unemp_levels), collapse = ", "),
    ". Keep the Wave-16 factor levels consistent before Bayesian export.",
    call. = FALSE
  )
}

dat <- dat %>%
  mutate(
    state_id = droplevels(factor(.data$state_id)),
    unemp_code = as.integer(factor(.data$unemployment_duration_cat, levels = unemp_levels)),
    y = as.integer(.data$curr_empl),
    id = as.integer(.data$state_id)
  )

core_no_missing <- c("y", "id")
if (anyNA(dat[, core_no_missing, drop = FALSE])) {
  stop("Outcome y or grouping id contains missing values. These are not imputed.", call. = FALSE)
}
if (!all(dat$y %in% c(0L, 1L))) {
  stop("Outcome y must be coded 0/1.", call. = FALSE)
}

numeric_terms <- c(
  "age_z", "age_z2", "female",
  "years_education", "ict_wle", "ict_meta", "math_wle", "reading_wle",
  "unemp_code",
  "migration_background", "father_unemployed", "mother_unemployed"
)

for (v in numeric_terms) {
  dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
}

bayes_input <- dat %>%
  transmute(
    row_id = row_number(),
    ID_t = as.numeric(.data$ID_t),
    y = as.integer(.data$y),
    id = as.integer(.data$id),
    age_z = .data$age_z,
    age_z2 = .data$age_z2,
    female = .data$female,
    years_education = .data$years_education,
    ict_wle = .data$ict_wle,
    ict_meta = .data$ict_meta,
    math_wle = .data$math_wle,
    reading_wle = .data$reading_wle,
    unemp_code = .data$unemp_code,
    migration_background = .data$migration_background,
    father_unemployed = .data$father_unemployed,
    mother_unemployed = .data$mother_unemployed
  )

bayes_input_file <- file.path(input_dir, "sc6_wave16_bayes_input.csv")
readr::write_csv(bayes_input, bayes_input_file, na = "")

state_lookup <- dat %>%
  distinct(id, state_id, FederalState) %>%
  arrange(.data$id) %>%
  transmute(id = .data$id, state_id = as.character(.data$state_id), FederalState = as.character(.data$FederalState))
readr::write_csv(state_lookup, file.path(input_dir, "sc6_wave16_state_lookup.csv"))

unemp_lookup <- tibble(
  unemp_code = seq_along(unemp_levels),
  unemployment_duration_cat = unemp_levels,
  reference_in_models = unemp_code == 1L
)
readr::write_csv(unemp_lookup, file.path(input_dir, "sc6_wave16_unemployment_lookup.csv"))

model_terms <- tibble::tribble(
  ~model, ~block, ~terms,
  "M1", "Baseline", "age_z,age_z2,female",
  "M2", "Baseline + human capital", "age_z,age_z2,female,years_education,ict_wle,ict_meta,math_wle,reading_wle",
  "M3", "Baseline + labour-market history", "age_z,age_z2,female,unemp_lt1,unemp_1to2,unemp_2to5,unemp_ge5",
  "M4", "Baseline + social background", "age_z,age_z2,female,migration_background,father_unemployed,mother_unemployed",
  "M5", "Baseline + all substantive blocks", "age_z,age_z2,female,years_education,ict_wle,ict_meta,math_wle,reading_wle,unemp_lt1,unemp_1to2,unemp_2to5,unemp_ge5,migration_background,father_unemployed,mother_unemployed"
)
readr::write_csv(model_terms, file.path(input_dir, "sc6_wave16_bayes_model_terms.csv"))

imputation_terms <- tibble(
  role = c(rep("baseline", 3), rep("human_capital", 5), "labour_history", rep("social_background", 3)),
  variable = numeric_terms
)
readr::write_csv(imputation_terms, file.path(input_dir, "sc6_wave16_bayes_imputation_terms.csv"))

audit <- tibble(
  variable = names(bayes_input),
  n = nrow(bayes_input),
  n_missing = vapply(bayes_input, function(x) sum(is.na(x)), integer(1)),
  n_observed = vapply(bayes_input, function(x) sum(!is.na(x)), integer(1)),
  n_unique_observed = vapply(bayes_input, function(x) length(unique(x[!is.na(x)])), integer(1)),
  observed_values = vapply(bayes_input, function(x) {
    vals <- unique(as.character(x[!is.na(x)]))
    paste(vals[seq_len(min(length(vals), 10L))], collapse = " | ")
  }, character(1))
) %>%
  mutate(
    imputed_in_julia = .data$variable %in% numeric_terms & .data$n_missing > 0L,
    problematic_for_cartbb = .data$imputed_in_julia & (.data$n_observed == 0L | .data$n_unique_observed < 2L)
  )
readr::write_csv(audit, file.path(input_dir, "sc6_wave16_bayes_input_audit.csv"))

bad <- audit %>% filter(.data$problematic_for_cartbb)
if (nrow(bad) > 0L) {
  stop(
    "At least one Bayesian imputation target is all-missing or constant: ",
    paste(bad$variable, collapse = ", "),
    ". Inspect sc6_wave16_bayes_input_audit.csv.",
    call. = FALSE
  )
}

# Compose Julia commands --------------------------------------------------

script_dir <- file.path(SC6_SCRIPT_DIR, "neps_sc6_empirical")
julia_script <- file.path(script_dir, "71_SC6_run_bayes_medoid.jl")
module_file <- file.path(SC6_SCRIPT_DIR, "julia", "src", "BayesREProbit.jl")

if (!file.exists(julia_script)) {
  warning("Julia runner not found yet: ", julia_script, call. = FALSE)
}
if (!file.exists(module_file)) {
  warning("Julia module not found yet: ", module_file, call. = FALSE)
}

ctl <- SC6_CONTROL
ctl$medoid_n_iter <- ctl$medoid_n_iter %||% ctl$bayes_n_iter %||% 20000L
ctl$medoid_burnin <- ctl$medoid_burnin %||% max(1000L, floor((ctl$medoid_n_iter %||% 20000L) / 4))
ctl$medoid_thin <- ctl$medoid_thin %||% 10L
ctl$bayes_n_iter <- ctl$bayes_n_iter %||% 20000L
ctl$bayes_burnin <- ctl$bayes_burnin %||% 4000L
ctl$bayes_thin <- ctl$bayes_thin %||% 2L
ctl$chib_reduced_iter <- ctl$chib_reduced_iter %||% 4000L
ctl$chib_reduced_burnin <- ctl$chib_reduced_burnin %||% 1000L
ctl$chib_reduced_thin <- ctl$chib_reduced_thin %||% 2L
ctl$cart_every <- ctl$cart_every %||% 10L
ctl$cart_max_depth <- ctl$cart_max_depth %||% 5L
ctl$cart_min_samples_leaf <- ctl$cart_min_samples_leaf %||% 10L
ctl$cart_max_donors <- ctl$cart_max_donors %||% 50L
ctl$tau2_beta <- ctl$tau2_beta %||% 1.0
ctl$c0 <- ctl$c0 %||% 2.0
ctl$d0 <- ctl$d0 %||% 1.0
ctl$seed <- ctl$seed %||% 20260628L

make_julia_args <- function(smoke = FALSE) {
  if (smoke) {
    medoid_n_iter <- 800L; medoid_burnin <- 200L; medoid_thin <- 10L
    n_iter <- 1200L; burnin <- 300L; thin <- 5L
    chib_iter <- 800L; chib_burnin <- 200L; chib_thin <- 2L
    fixed_points <- "mean"
  } else {
    medoid_n_iter <- ctl$medoid_n_iter; medoid_burnin <- ctl$medoid_burnin; medoid_thin <- ctl$medoid_thin
    n_iter <- ctl$bayes_n_iter; burnin <- ctl$bayes_burnin; thin <- ctl$bayes_thin
    chib_iter <- ctl$chib_reduced_iter; chib_burnin <- ctl$chib_reduced_burnin; chib_thin <- ctl$chib_reduced_thin
    fixed_points <- Sys.getenv("SC6_BAYES_FIXED_POINTS", unset = "mean")
  }

  c(
    "--data", normalizePath(bayes_input_file, winslash = "/", mustWork = FALSE),
    "--outdir", normalizePath(SC6_BAYES_DIR, winslash = "/", mustWork = FALSE),
    "--module", normalizePath(module_file, winslash = "/", mustWork = FALSE),
    "--medoid_reference_model", ctl$medoid_reference_model %||% "M5",
    "--medoid_n_iter", as.character(medoid_n_iter),
    "--medoid_burnin", as.character(medoid_burnin),
    "--medoid_thin", as.character(medoid_thin),
    "--n_iter", as.character(n_iter),
    "--burnin", as.character(burnin),
    "--thin", as.character(thin),
    "--chib_reduced_iter", as.character(chib_iter),
    "--chib_reduced_burnin", as.character(chib_burnin),
    "--chib_reduced_thin", as.character(chib_thin),
    "--cart_every", as.character(ctl$cart_every),
    "--cart_max_depth", as.character(ctl$cart_max_depth),
    "--cart_min_samples_leaf", as.character(ctl$cart_min_samples_leaf),
    "--cart_max_donors", as.character(ctl$cart_max_donors),
    "--tau2_beta", as.character(ctl$tau2_beta),
    "--c0", as.character(ctl$c0),
    "--d0", as.character(ctl$d0),
    "--seed", as.character(ctl$seed),
    "--fixed_points", fixed_points,
    "--verbose_every", Sys.getenv("SC6_BAYES_VERBOSE_EVERY", unset = "1000")
  )
}

write_command_files <- function(args, suffix) {
  qargs <- vapply(args, shQuote, character(1))
  qscript <- shQuote(normalizePath(julia_script, winslash = "/", mustWork = FALSE))
  cmd <- paste(c("julia", "--project=.", qscript, qargs), collapse = " ")
  writeLines(cmd, file.path(input_dir, paste0("run_sc6_bayes_", suffix, ".sh")))
  writeLines(cmd, file.path(input_dir, paste0("run_sc6_bayes_", suffix, ".bat")))
  invisible(cmd)
}

smoke_args <- make_julia_args(smoke = TRUE)
paper_args <- make_julia_args(smoke = FALSE)
smoke_cmd <- write_command_files(smoke_args, "smoke")
paper_cmd <- write_command_files(paper_args, "paper")

run_julia <- tolower(Sys.getenv("SC6_RUN_JULIA", unset = "false")) %in% c("true", "1", "yes")
julia_executed <- FALSE
julia_status <- "not requested"

if (run_julia) {
  if (!file.exists(julia_script)) stop("Cannot run Julia; runner not found: ", julia_script, call. = FALSE)
  if (!file.exists(module_file)) stop("Cannot run Julia; module not found: ", module_file, call. = FALSE)
  julia_bin <- Sys.getenv("JULIA_BIN", unset = "julia")
  message_ts("Launching Julia Bayesian CART-BB workflow.")
  status <- system2(julia_bin, args = c("--project=.", normalizePath(julia_script, winslash = "/", mustWork = FALSE), paper_args),
                    stdout = file.path(SC6_BAYES_DIR, "julia_stdout.log"),
                    stderr = file.path(SC6_BAYES_DIR, "julia_stderr.log"))
  julia_executed <- TRUE
  julia_status <- paste0("system2 status=", status)
  if (!identical(status, 0L)) {
    stop("Julia Bayesian workflow returned non-zero status. See julia_stdout.log and julia_stderr.log in ", SC6_BAYES_DIR, call. = FALSE)
  }
}

status_tbl <- tibble(
  input_csv = normalizePath(bayes_input_file, winslash = "/", mustWork = FALSE),
  bayes_dir = normalizePath(SC6_BAYES_DIR, winslash = "/", mustWork = FALSE),
  julia_script = normalizePath(julia_script, winslash = "/", mustWork = FALSE),
  module_file = normalizePath(module_file, winslash = "/", mustWork = FALSE),
  smoke_command_file = normalizePath(file.path(input_dir, "run_sc6_bayes_smoke.sh"), winslash = "/", mustWork = FALSE),
  paper_command_file = normalizePath(file.path(input_dir, "run_sc6_bayes_paper.sh"), winslash = "/", mustWork = FALSE),
  julia_executed = julia_executed,
  julia_status = julia_status
)
readr::write_csv(status_tbl, file.path(input_dir, "sc6_bayes_export_status.csv"))

message_ts("Bayesian Julia input written to: ", bayes_input_file)
message_ts("Julia smoke command: ", smoke_cmd)
message_ts("Julia paper command: ", paper_cmd)
message_ts("Set Sys.setenv(SC6_RUN_JULIA='true') if you want this R stage to launch Julia automatically.")
