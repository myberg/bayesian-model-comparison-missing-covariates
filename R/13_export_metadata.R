find_index_file_column <- function(idx) {
  candidates <- c("file_mis", "file_missing", "missing_file", "file")
  out <- intersect(candidates, names(idx))
  if (length(out) == 0L) {
    stop("No missing-data file column found in simulation index.")
  }
  out[1]
}

export_missingness_metadata <- function(index_file,
                                        out_dir,
                                        covariates = c("x2", "x3", "x5B", "x5C", "x7", "x10")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  idx <- read.csv(index_file)
  file_col <- find_index_file_column(idx)
  
  miss_long <- lapply(seq_len(nrow(idx)), function(i) {
    dat <- read.csv(idx[[file_col]][i])
    vars <- intersect(covariates, names(dat))
    
    data.frame(
      rep = idx$rep_id[i],
      mechanism = idx$mechanism[i],
      variable = vars,
      missing_rate = vapply(vars, function(v) mean(is.na(dat[[v]])), numeric(1)),
      stringsAsFactors = FALSE
    )
  })
  
  miss_long <- dplyr::bind_rows(miss_long)
  
  miss_summary <- miss_long |>
    dplyr::group_by(.data$mechanism, .data$variable) |>
    dplyr::summarise(
      mean_missing_rate = mean(.data$missing_rate, na.rm = TRUE),
      sd_missing_rate = stats::sd(.data$missing_rate, na.rm = TRUE),
      min_missing_rate = min(.data$missing_rate, na.rm = TRUE),
      max_missing_rate = max(.data$missing_rate, na.rm = TRUE),
      n_rep = dplyr::n_distinct(.data$rep),
      .groups = "drop"
    )
  
  readr::write_csv(
    miss_long,
    file.path(out_dir, "realized_missingness_rates_by_replication.csv")
  )
  
  readr::write_csv(
    miss_summary,
    file.path(out_dir, "realized_missingness_rates.csv")
  )
  
  invisible(miss_summary)
}

export_complete_case_rates <- function(index_file,
                                       out_dir,
                                       complete_case_vars = c("y", "x2", "x3", "x4",
                                                              "x5B", "x5C",
                                                              "x7", "x8B", "x8C", "x10")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  idx <- read.csv(index_file)
  file_col <- find_index_file_column(idx)
  
  cc <- lapply(seq_len(nrow(idx)), function(i) {
    dat <- read.csv(idx[[file_col]][i])
    vars <- intersect(complete_case_vars, names(dat))
    
    cc_vec <- stats::complete.cases(dat[, vars, drop = FALSE])
    
    data.frame(
      rep = idx$rep_id[i],
      mechanism = idx$mechanism[i],
      n_total = nrow(dat),
      n_complete = sum(cc_vec),
      complete_case_rate = mean(cc_vec),
      stringsAsFactors = FALSE
    )
  })
  
  cc <- dplyr::bind_rows(cc)
  
  cc_summary <- cc |>
    dplyr::group_by(.data$mechanism) |>
    dplyr::summarise(
      mean_complete_case_rate = mean(.data$complete_case_rate, na.rm = TRUE),
      sd_complete_case_rate = stats::sd(.data$complete_case_rate, na.rm = TRUE),
      mean_n_complete = mean(.data$n_complete, na.rm = TRUE),
      n_rep = dplyr::n_distinct(.data$rep),
      .groups = "drop"
    )
  
  readr::write_csv(
    cc,
    file.path(out_dir, "complete_case_rates_by_replication.csv")
  )
  
  readr::write_csv(
    cc_summary,
    file.path(out_dir, "complete_case_rates.csv")
  )
  
  invisible(cc_summary)
}