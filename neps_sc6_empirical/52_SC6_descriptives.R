###
# 52_SC6_descriptives.R
#
# Descriptive summaries for the prepared NEPS SC6 empirical analysis dataset.
# This script is called by 00_master_SC6_empirical_analysis.R.
###

library(dplyr)

processed_dir <- Sys.getenv("SC6_OUTPUT_DIR", unset = file.path(getwd(), "data", "processed"))
desc_dir <- Sys.getenv("SC6_DESC_DIR", unset = file.path(getwd(), "results", "sc6_empirical", "descriptives"))

dir.create(desc_dir, recursive = TRUE, showWarnings = FALSE)

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

sample_flow_file <- file.path(processed_dir, "SC6_wave16_sample_flow.csv")
if (file.exists(sample_flow_file)) {
  sample_flow <- read.csv(sample_flow_file)
  utils::write.csv(sample_flow, file.path(desc_dir, "sample_flow.csv"), row.names = FALSE)
} else {
  warning("Sample-flow file not found: ", sample_flow_file, call. = FALSE)
}

write_csv <- function(x, file) {
  utils::write.csv(x, file.path(desc_dir, file), row.names = FALSE)
}

continuous_vars <- intersect(
  c(
    "age", "age_z", "age_z2", "years_education",
    "ict_wle", "ict_meta", "math_wle", "reading_wle",
    "unemployment_months", "unemployment_years"
  ),
  names(dat)
)

categorical_vars <- intersect(
  c(
    "curr_empl", "female", "FederalState", "CASMIN",
    "unemployment_duration_cat", "migration_background",
    "father_unemployed", "mother_unemployed", "university_degree"
  ),
  names(dat)
)

continuous_summary <- lapply(continuous_vars, function(v) {
  x <- dat[[v]]
  tibble::tibble(
    variable = v,
    n_observed = sum(!is.na(x)),
    n_missing = sum(is.na(x)),
    missing_rate = mean(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    q25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)),
    median = stats::median(x, na.rm = TRUE),
    q75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}) %>%
  bind_rows()

categorical_summary <- lapply(categorical_vars, function(v) {
  dat %>%
    count(level = as.character(.data[[v]]), name = "n") %>%
    mutate(
      variable = v,
      level = if_else(is.na(.data$level), "<missing>", .data$level),
      n_total = sum(.data$n),
      n_observed = sum(.data$n[.data$level != "<missing>"]),
      proportion_total = .data$n / .data$n_total,
      proportion_observed = if_else(
        .data$level == "<missing>" | .data$n_observed == 0L,
        NA_real_,
        .data$n / .data$n_observed
      )
    ) %>%
    select(variable, level, n, n_total, n_observed,
           proportion_total, proportion_observed)
}) %>%
  bind_rows()

missingness_summary <- dat %>%
  summarise(across(everything(), ~ mean(is.na(.x)))) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  rename(missing_rate = V1) %>%
  arrange(desc(.data$missing_rate), .data$variable)

model_spec_table <- lapply(names(model_variables), function(model_name) {
  tibble::tibble(
    model = model_name,
    block = case_when(
      model_name == "M1" ~ "Baseline",
      model_name == "M2" ~ "Baseline + human capital",
      model_name == "M3" ~ "Baseline + labour-market history",
      model_name == "M4" ~ "Baseline + social background",
      model_name == "M5" ~ "Baseline + all substantive blocks",
      TRUE ~ model_name
    ),
    fixed_effect_variables = paste(model_variables[[model_name]], collapse = ", "),
    random_effect = "Federal-state random intercept"
  )
}) %>%
  bind_rows()

complete_case_by_model <- lapply(names(model_variables), function(model_name) {
  vars <- unique(c("curr_empl", "state_id", model_variables[[model_name]]))
  n_complete <- sum(stats::complete.cases(dat[, vars, drop = FALSE]))
  n_total <- nrow(dat)
  tibble::tibble(
    model = model_name,
    n_complete = n_complete,
    n_total = n_total,
    complete_rate = n_complete / n_total
  )
}) %>%
  bind_rows()

write_csv(continuous_summary, "continuous_summary.csv")
write_csv(categorical_summary, "categorical_summary.csv")
write_csv(missingness_summary, "missingness_summary.csv")
write_csv(model_spec_table, "model_specification_table.csv")
write_csv(complete_case_by_model, "complete_case_by_model.csv")

message("Descriptive summaries written to: ", normalizePath(desc_dir, mustWork = FALSE))
