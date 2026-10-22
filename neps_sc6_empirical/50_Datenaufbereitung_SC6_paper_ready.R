###
# NEPS SC6 empirical illustration
# Data preparation for random-effects probit model comparison
#
# Purpose:
#   Build one analytic Wave-16 employment dataset that is consistent with the
#   paper design: fixed outcome and grouping variables, incomplete substantive
#   covariates retained for imputation/model comparison.
#
# Required input:
#   NEPS SC6 SUF 16.0.0 Stata files in one directory.
#
# Configure paths either here or via environment variables:
#   NEPS_SC6_STATA_DIR=/path/to/Stata14
#   SC6_OUTPUT_DIR=/path/to/processed/output
###

library(dplyr)
library(haven)
library(lubridate)

# Paths -------------------------------------------------------------------

neps_dir <- Sys.getenv(
  "NEPS_SC6_STATA_DIR",
  unset = file.path(getwd(), "data", "neps_sc6", "raw", "Stata14")
)

out_dir <- Sys.getenv(
  "SC6_OUTPUT_DIR",
  unset = file.path(getwd(), "data", "neps_sc6", "processed")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_neps <- function(file) {
  path <- file.path(neps_dir, file)
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  haven::read_stata(path)
}

missing_codes <- c(-98, -97, -95, -94, -90, -56, -55, -54, -20)

clean_num <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  x_num[x_num %in% missing_codes] <- NA_real_
  x_num
}

safe_left_join <- function(x, y, by = "ID_t") {
  stopifnot(!anyDuplicated(y[[by]]))
  dplyr::left_join(x, y, by = by)
}

wave16_or_latest <- function(dat, wave_value = 16) {
  if ("wave" %in% names(dat)) {
    out <- dat %>%
      filter(.data$wave == wave_value) %>%
      arrange(.data$ID_t) %>%
      distinct(.data$ID_t, .keep_all = TRUE)
    
    if (nrow(out) > 0) {
      return(out)
    }
    
    return(dat %>%
             arrange(.data$ID_t, desc(.data$wave)) %>%
             distinct(.data$ID_t, .keep_all = TRUE))
  }
  
  dat %>%
    arrange(.data$ID_t) %>%
    distinct(.data$ID_t, .keep_all = TRUE)
}

make_binary <- function(x, one_values, zero_values) {
  dplyr::case_when(
    x %in% one_values ~ 1L,
    x %in% zero_values ~ 0L,
    TRUE ~ NA_integer_
  )
}

# Load source files -------------------------------------------------------

ptarget_raw <- read_neps("SC6_pTarget_D_16-0-0.dta")
basics_raw  <- read_neps("SC6_Basics_D_16-0-0.dta")
comp_raw    <- read_neps("SC6_xTargetCompetencies_D_16-0-0.dta")
unemp_raw   <- read_neps("SC6_spUnemp_D_16-0-0.dta")

ptarget_w16 <- ptarget_raw %>%
  filter(.data$wave == 16) %>%
  arrange(.data$ID_t) %>%
  distinct(.data$ID_t, .keep_all = TRUE)

wave16_ids <- ptarget_w16 %>%
  transmute(ID_t = as.numeric(.data$ID_t))

basics_w16 <- basics_raw %>%
  semi_join(wave16_ids, by = "ID_t") %>%
  wave16_or_latest(wave_value = 16)

comp_wide <- comp_raw %>%
  semi_join(wave16_ids, by = "ID_t") %>%
  wave16_or_latest(wave_value = 16)

# Outcome and baseline variables ----------------------------------------

outcome_base <- basics_w16 %>%
  transmute(
    ID_t = as.numeric(.data$ID_t),
    curr_empl_raw = clean_num(.data$tx29060),
    curr_empl = make_binary(.data$curr_empl_raw, one_values = c(1), zero_values = c(0, 2)),
    curr_unempl_raw = clean_num(.data$tx29080),
    federal_state_code = clean_num(.data$t751001_g2),
    FederalState = dplyr::recode(
      as.character(.data$federal_state_code),
      `1`  = "Schleswig-Holstein",
      `2`  = "Hamburg",
      `3`  = "Lower Saxony",
      `4`  = "Bremen",
      `5`  = "North Rhine-Westphalia",
      `6`  = "Hesse",
      `7`  = "Rhineland-Palatinate",
      `8`  = "Baden-Wuerttemberg",
      `9`  = "Bavaria",
      `10` = "Saarland",
      `11` = "Berlin",
      `12` = "Brandenburg",
      `13` = "Mecklenburg-Western Pomerania",
      `14` = "Saxony",
      `15` = "Saxony-Anhalt",
      `16` = "Thuringia",
      .default = NA_character_
    ),
    state_id = factor(.data$FederalState)
  )

# Age and gender are Wave-16 variables. The parental background indicators are
# time-stable/background variables and may not be observed in Wave 16 for all
# respondents. We therefore retrieve the most recent non-missing value across all
# available waves for Wave-16 respondents. This avoids turning a valid background
# variable into an all-missing or constant Wave-16 artifact.
latest_nonmissing_value <- function(dat, source_var, out_name) {
  dat %>%
    semi_join(wave16_ids, by = "ID_t") %>%
    transmute(
      ID_t = as.numeric(.data$ID_t),
      wave = if ("wave" %in% names(dat)) suppressWarnings(as.numeric(.data$wave)) else 0,
      value = clean_num(.data[[source_var]])
    ) %>%
    arrange(.data$ID_t, desc(.data$wave)) %>%
    group_by(.data$ID_t) %>%
    summarise(
      "{out_name}" := {
        vv <- .data$value[!is.na(.data$value)]
        if (length(vv) == 0L) NA_real_ else vv[1]
      },
      .groups = "drop"
    )
}

parent_background <- latest_nonmissing_value(ptarget_raw, "t731451", "father_unemployed_raw") %>%
  left_join(
    latest_nonmissing_value(ptarget_raw, "t731401", "mother_unemployed_raw"),
    by = "ID_t"
  ) %>%
  mutate(
    father_unemployed = make_binary(.data$father_unemployed_raw, one_values = c(2), zero_values = c(1)),
    mother_unemployed = make_binary(.data$mother_unemployed_raw, one_values = c(2), zero_values = c(1))
  )

personal_base <- ptarget_w16 %>%
  transmute(
    ID_t = as.numeric(.data$ID_t),
    birth_year = clean_num(.data$t70000y),
    age = 2022 - .data$birth_year,
    female = make_binary(clean_num(.data$t700001), one_values = c(2), zero_values = c(1))
  ) %>%
  left_join(parent_background, by = "ID_t")

# Human capital variables -------------------------------------------------

human_capital <- basics_w16 %>%
  transmute(
    ID_t = as.numeric(.data$ID_t),
    casmin_raw = clean_num(.data$tx28101),
    CASMIN = dplyr::case_when(
      is.na(.data$casmin_raw) ~ NA_character_,
      .data$casmin_raw %in% c(0, 1, 3) ~ "1a/1b/2b",
      .data$casmin_raw %in% c(2, 4) ~ "1c/2a",
      .data$casmin_raw %in% c(5, 6) ~ "2c",
      .data$casmin_raw %in% c(7, 8) ~ "3a/3b",
      TRUE ~ NA_character_
    ),
    CASMIN = factor(.data$CASMIN, levels = c("1a/1b/2b", "1c/2a", "2c", "3a/3b")),
    years_education = clean_num(.data$tx28102),
    isced_raw = clean_num(.data$tx28103),
    university_degree = dplyr::case_when(
      is.na(.data$isced_raw) ~ NA_integer_,
      .data$isced_raw %in% c(8, 9, 10) ~ 1L,
      TRUE ~ 0L
    )
  )

competencies <- comp_wide %>%
  transmute(
    ID_t = as.numeric(.data$ID_t),
    ict_wle = clean_num(.data$ica5_sc1),
    ict_meta = clean_num(.data$mpa5ic_sc6),
    math_wle = clean_num(.data$maa9_sc1),
    reading_wle = clean_num(.data$rea9_sc1)
  )

# Labour-market history ---------------------------------------------------

# The spell file contains unemployment episodes. Respondents without a record
# are treated as having no previous unemployment spell. Incomplete spells remain
# missing rather than being turned into zero duration.
analysis_end_date <- as.Date("2023-05-01")

unemp_spells <- unemp_raw %>%
  semi_join(wave16_ids, by = "ID_t") %>%
  mutate(
    ID_t = as.numeric(.data$ID_t),
    start_year = clean_num(.data$ts2511y),
    start_month = clean_num(.data$ts2511m),
    end_year = clean_num(.data$ts2512y),
    end_month = clean_num(.data$ts2512m),
    start_date = lubridate::make_date(.data$start_year, .data$start_month, 1),
    end_date = lubridate::make_date(.data$end_year, .data$end_month, 1),
    end_date = if_else(is.na(.data$end_date), analysis_end_date, .data$end_date),
    unemp_months_spell = interval(.data$start_date, .data$end_date) %/% months(1),
    unemp_months_spell = if_else(.data$unemp_months_spell < 0, NA_integer_, .data$unemp_months_spell)
  )

unemployment_history <- unemp_spells %>%
  group_by(.data$ID_t) %>%
  summarise(
    unemp_spell_count = n(),
    any_incomplete_unemp_spell = any(is.na(.data$unemp_months_spell)),
    unemployment_months = if_else(
      any(is.na(.data$unemp_months_spell)),
      NA_real_,
      as.numeric(sum(.data$unemp_months_spell, na.rm = TRUE))
    ),
    .groups = "drop"
  ) %>%
  right_join(wave16_ids, by = "ID_t") %>%
  mutate(
    unemp_spell_count = if_else(is.na(.data$unemp_spell_count), 0L, .data$unemp_spell_count),
    any_incomplete_unemp_spell = if_else(
      is.na(.data$any_incomplete_unemp_spell),
      FALSE,
      .data$any_incomplete_unemp_spell
    ),
    unemployment_months = if_else(.data$unemp_spell_count == 0L, 0, .data$unemployment_months),
    unemployment_years = .data$unemployment_months / 12,
    unemployment_duration_cat = case_when(
      is.na(.data$unemployment_months) ~ NA_character_,
      .data$unemployment_months == 0 ~ "no previous spell",
      .data$unemployment_months < 12 ~ "<1 year",
      .data$unemployment_months < 24 ~ "1 to <2 years",
      .data$unemployment_months < 60 ~ "2 to <5 years",
      .data$unemployment_months >= 60 ~ ">=5 years"
    ),
    unemployment_duration_cat = factor(
      .data$unemployment_duration_cat,
      levels = c("no previous spell", "<1 year", "1 to <2 years", "2 to <5 years", ">=5 years")
    )
  )

# Social-background variables -------------------------------------------

social_background <- basics_w16 %>%
  transmute(
    ID_t = as.numeric(.data$ID_t),
    migration_generation_raw = clean_num(.data$t400500_g1),
    migration_background = case_when(
      is.na(.data$migration_generation_raw) ~ NA_integer_,
      .data$migration_generation_raw %in% 1:6 ~ 1L,
      TRUE ~ 0L
    )
  )

# Assemble analytic data --------------------------------------------------

sc6_analysis_all <- wave16_ids %>%
  safe_left_join(outcome_base, by = "ID_t") %>%
  safe_left_join(personal_base, by = "ID_t") %>%
  safe_left_join(human_capital, by = "ID_t") %>%
  safe_left_join(competencies, by = "ID_t") %>%
  safe_left_join(unemployment_history, by = "ID_t") %>%
  safe_left_join(social_background, by = "ID_t") %>%
  mutate(
    state_id = droplevels(factor(.data$FederalState))
  )

# Sequential sample-flow table.
# Exclusions are counted step by step, so a respondent removed earlier is not
# counted again in a later exclusion category.
flow0 <- sc6_analysis_all
flow1 <- flow0 %>% filter(!is.na(.data$age))
flow2 <- flow1 %>% filter(.data$age < 65)
flow3 <- flow2 %>% filter(!is.na(.data$curr_empl))
flow4 <- flow3 %>% filter(!is.na(.data$FederalState))

sample_flow <- tibble::tibble(
  step = c(
    "Wave 16 respondents",
    "Missing interview age",
    "Age 65 or older at interview",
    "Missing employment outcome",
    "Missing federal state of residence",
    "Analytic sample"
  ),
  excluded = c(
    NA_integer_,
    nrow(flow0) - nrow(flow1),
    nrow(flow1) - nrow(flow2),
    nrow(flow2) - nrow(flow3),
    nrow(flow3) - nrow(flow4),
    NA_integer_
  ),
  remaining = c(
    nrow(flow0),
    nrow(flow1),
    nrow(flow2),
    nrow(flow3),
    nrow(flow4),
    nrow(flow4)
  )
)

# Non-exclusive diagnostics are useful for checks only; do not report them as
# the main sample-flow table because exclusion reasons may overlap.
sample_flow_diagnostics_nonexclusive <- tibble::tibble(
  criterion = c(
    "Wave 16 respondents",
    "Missing interview age",
    "Age 65 or older at interview",
    "Missing employment outcome",
    "Missing federal state of residence"
  ),
  n = c(
    nrow(flow0),
    sum(is.na(flow0$age)),
    sum(!is.na(flow0$age) & flow0$age >= 65),
    sum(is.na(flow0$curr_empl)),
    sum(is.na(flow0$FederalState))
  ),
  proportion = .data$n / nrow(flow0)
)


sc6_analysis <- flow4 %>%
  mutate(
    state_id = droplevels(.data$state_id),
    age_z = as.numeric(scale(.data$age)),
    age_z2_raw = .data$age_z^2,
    age_z2 = as.numeric(scale(.data$age_z2_raw))
  )

missingness <- sc6_analysis %>%
  summarise(across(everything(), ~ mean(is.na(.x)))) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "missing_rate") %>%
  arrange(desc(.data$missing_rate), .data$variable)

state_counts <- sc6_analysis %>%
  count(.data$FederalState, .data$state_id, name = "n") %>%
  arrange(.data$state_id)

# Model-variable sets -----------------------------------------------------

baseline_vars <- c("age_z", "age_z2", "female")
human_capital_vars <- c("years_education", "ict_wle", "ict_meta", "math_wle", "reading_wle")
labour_history_vars <- c("unemployment_duration_cat")
social_background_vars <- c("migration_background", "father_unemployed", "mother_unemployed")

model_variables <- list(
  M1 = baseline_vars,
  M2 = c(baseline_vars, human_capital_vars),
  M3 = c(baseline_vars, labour_history_vars),
  M4 = c(baseline_vars, social_background_vars),
  M5 = c(baseline_vars, human_capital_vars, labour_history_vars, social_background_vars)
)

model_formulas <- list(
  M1 = curr_empl ~ age_z + age_z2 + female + (1 | state_id),
  M2 = curr_empl ~ age_z + age_z2 + female +
    years_education + ict_wle + ict_meta + math_wle + reading_wle + (1 | state_id),
  M3 = curr_empl ~ age_z + age_z2 + female +
    unemployment_duration_cat + (1 | state_id),
  M4 = curr_empl ~ age_z + age_z2 + female +
    migration_background + father_unemployed + mother_unemployed + (1 | state_id),
  M5 = curr_empl ~ age_z + age_z2 + female +
    years_education + ict_wle + ict_meta + math_wle + reading_wle +
    unemployment_duration_cat +
    migration_background + father_unemployed + mother_unemployed + (1 | state_id)
)

# Julia/Gibbs helper objects ---------------------------------------------

state_id_numeric <- as.integer(sc6_analysis$state_id)
D_N <- model.matrix(~ 0 + sc6_analysis$state_id)
colnames(D_N) <- paste0("state_", seq_len(ncol(D_N)))

julia_data <- sc6_analysis %>%
  mutate(
    Intercept = 1,
    state_id_numeric = state_id_numeric,
    unemployment_duration_num = as.integer(.data$unemployment_duration_cat),
    CASMIN_num = as.integer(.data$CASMIN)
  )

# Save outputs ------------------------------------------------------------

saveRDS(sc6_analysis, file.path(out_dir, "SC6_wave16_employment_analysis.rds"))
saveRDS(julia_data, file.path(out_dir, "SC6_wave16_employment_analysis_julia.rds"))
saveRDS(D_N, file.path(out_dir, "SC6_wave16_state_design_matrix.rds"))
saveRDS(model_variables, file.path(out_dir, "SC6_wave16_model_variables.rds"))
saveRDS(model_formulas, file.path(out_dir, "SC6_wave16_model_formulas.rds"))

haven::write_dta(sc6_analysis, file.path(out_dir, "SC6_wave16_employment_analysis.dta"))
write.csv(sample_flow, file.path(out_dir, "SC6_wave16_sample_flow.csv"), row.names = FALSE)
write.csv(sample_flow_diagnostics_nonexclusive,
          file.path(out_dir, "SC6_wave16_sample_flow_diagnostics_nonexclusive.csv"),
          row.names = FALSE)
write.csv(missingness, file.path(out_dir, "SC6_wave16_missingness.csv"), row.names = FALSE)
write.csv(state_counts, file.path(out_dir, "SC6_wave16_state_counts.csv"), row.names = FALSE)

# Optional: also copy paper-relevant preparation outputs into the current run's
# descriptives folder if the master workflow has set SC6_DESC_DIR.
desc_dir <- Sys.getenv("SC6_DESC_DIR", unset = NA_character_)
if (!is.na(desc_dir) && nzchar(desc_dir)) {
  dir.create(desc_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(sample_flow, file.path(desc_dir, "sample_flow.csv"), row.names = FALSE)
  write.csv(sample_flow_diagnostics_nonexclusive,
            file.path(desc_dir, "sample_flow_diagnostics_nonexclusive.csv"),
            row.names = FALSE)
  write.csv(missingness, file.path(desc_dir, "missingness_from_preparation.csv"), row.names = FALSE)
  write.csv(state_counts, file.path(desc_dir, "state_counts.csv"), row.names = FALSE)
}

message("Analytic sample size: ", nrow(sc6_analysis))
message("Number of federal states: ", nlevels(sc6_analysis$state_id))
message("Outputs written to: ", normalizePath(out_dir, mustWork = FALSE))
