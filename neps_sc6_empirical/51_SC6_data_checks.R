###
# 51_SC6_data_checks.R
#
# Checks for the prepared NEPS SC6 empirical analysis dataset.
# This script is called by 00_master_SC6_empirical_analysis.R.
###

library(dplyr)

processed_dir <- Sys.getenv("SC6_OUTPUT_DIR", unset = file.path(getwd(), "data", "processed"))
check_dir <- Sys.getenv("SC6_CHECK_DIR", unset = file.path(getwd(), "results", "sc6_empirical", "checks"))

dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

analysis_file <- file.path(processed_dir, "SC6_wave16_employment_analysis.rds")
model_variables_file <- file.path(processed_dir, "SC6_wave16_model_variables.rds")

if (!file.exists(analysis_file)) {
  stop("Prepared analysis file not found: ", analysis_file, call. = FALSE)
}
if (!file.exists(model_variables_file)) {
  stop("Model variable file not found: ", model_variables_file, call. = FALSE)
}

dat <- readRDS(analysis_file)
model_variables <- readRDS(model_variables_file)

write_csv <- function(x, file) {
  utils::write.csv(x, file.path(check_dir, file), row.names = FALSE)
}

required_core <- c("ID_t", "curr_empl", "FederalState", "state_id", "age", "age_z", "age_z2")
missing_core <- setdiff(required_core, names(dat))
if (length(missing_core) > 0L) {
  stop("Missing core variable(s): ", paste(missing_core, collapse = ", "), call. = FALSE)
}

all_model_vars <- unique(unlist(model_variables, use.names = FALSE))
missing_model_vars <- setdiff(all_model_vars, names(dat))
if (length(missing_model_vars) > 0L) {
  stop("Missing model variable(s): ", paste(missing_model_vars, collapse = ", "), call. = FALSE)
}

sample_check <- tibble::tibble(
  check = c(
    "n_rows",
    "n_unique_ID_t",
    "missing_outcome",
    "missing_state",
    "n_states",
    "min_state_n",
    "max_state_n",
    "min_age",
    "max_age"
  ),
  value = c(
    nrow(dat),
    dplyr::n_distinct(dat$ID_t),
    sum(is.na(dat$curr_empl)),
    sum(is.na(dat$state_id)),
    nlevels(droplevels(dat$state_id)),
    min(table(dat$state_id)),
    max(table(dat$state_id)),
    min(dat$age, na.rm = TRUE),
    max(dat$age, na.rm = TRUE)
  )
)

outcome_check <- dat %>%
  count(.data$curr_empl, name = "n") %>%
  mutate(proportion = .data$n / sum(.data$n))

state_check <- dat %>%
  count(.data$FederalState, .data$state_id, name = "n") %>%
  arrange(.data$state_id)

missingness_check <- dat %>%
  summarise(across(everything(), ~ mean(is.na(.x)))) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  rename(missing_rate = V1) %>%
  arrange(desc(.data$missing_rate), .data$variable)

complete_case_by_model <- lapply(names(model_variables), function(model_name) {
  vars <- unique(c("curr_empl", "state_id", model_variables[[model_name]]))
  n_complete <- sum(stats::complete.cases(dat[, vars, drop = FALSE]))
  n_total <- nrow(dat)
  tibble::tibble(
    model = model_name,
    n_variables = length(vars),
    n_complete = n_complete,
    n_total = n_total,
    complete_rate = n_complete / n_total
  )
}) %>%
  bind_rows()

variable_type_check <- tibble::tibble(
  variable = names(dat),
  class = vapply(dat, function(x) paste(class(x), collapse = "/"), character(1)),
  n_unique = vapply(dat, function(x) dplyr::n_distinct(x, na.rm = TRUE), integer(1)),
  missing_rate = vapply(dat, function(x) mean(is.na(x)), numeric(1))
)

range_check <- dat %>%
  summarise(
    age_min = min(.data$age, na.rm = TRUE),
    age_max = max(.data$age, na.rm = TRUE),
    ict_meta_min = if ("ict_meta" %in% names(dat)) min(.data$ict_meta, na.rm = TRUE) else NA_real_,
    ict_meta_max = if ("ict_meta" %in% names(dat)) max(.data$ict_meta, na.rm = TRUE) else NA_real_,
    unemployment_years_min = if ("unemployment_years" %in% names(dat)) min(.data$unemployment_years, na.rm = TRUE) else NA_real_,
    unemployment_years_max = if ("unemployment_years" %in% names(dat)) max(.data$unemployment_years, na.rm = TRUE) else NA_real_
  )

warnings <- character(0)

if (any(duplicated(dat$ID_t))) {
  warnings <- c(warnings, "ID_t is not unique.")
}
if (any(is.na(dat$curr_empl))) {
  warnings <- c(warnings, "Employment outcome contains missing values after analytic filtering.")
}
if (any(is.na(dat$state_id))) {
  warnings <- c(warnings, "State identifier contains missing values after analytic filtering.")
}
if (nlevels(droplevels(dat$state_id)) < 5L) {
  warnings <- c(warnings, "The grouping variable has few levels. Reconsider random-effects specification.")
}
if ("ict_meta" %in% names(dat)) {
  rng <- range(dat$ict_meta, na.rm = TRUE)
  if (is.finite(rng[1]) && is.finite(rng[2]) && (rng[1] < 0 || rng[2] > 1)) {
    warnings <- c(
      warnings,
      "ict_meta is labelled as a proportion but has values outside [0, 1]. Check whether it is standardized."
    )
  }
}

write_csv(sample_check, "sample_check.csv")
write_csv(outcome_check, "outcome_check.csv")
write_csv(state_check, "state_check.csv")
write_csv(missingness_check, "missingness_check.csv")
write_csv(complete_case_by_model, "complete_case_by_model.csv")
write_csv(variable_type_check, "variable_type_check.csv")
write_csv(range_check, "range_check.csv")

writeLines(
  if (length(warnings) == 0L) "No data-check warnings." else warnings,
  con = file.path(check_dir, "data_check_warnings.txt")
)

message("Data checks written to: ", normalizePath(check_dir, mustWork = FALSE))

