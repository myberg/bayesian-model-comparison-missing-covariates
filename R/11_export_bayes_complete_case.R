# ============================================================
# Export common complete-case datasets for Bayesian Julia runs
# ------------------------------------------------------------
# Adds file_cc and n_cc to data/simulated/simulation_index.csv.
# The common CC sample is defined over y, id, and all variables used
# in any candidate model so that Bayesian model comparisons use the
# same observations across all candidate models.
# ============================================================

.model_vars_from_specs_cc <- function(model_specs_arg) {
  unique(unlist(lapply(model_specs_arg, function(x) {
    all.vars(x$fixed)
  }), use.names = FALSE))
}

.resolve_index_path_cc <- function(p, index_file) {
  if (is.na(p) || !nzchar(p)) return(NA_character_)
  if (file.exists(p)) return(p)
  p2 <- file.path(dirname(index_file), basename(p))
  if (file.exists(p2)) return(p2)
  p
}

export_common_complete_case_for_julia <- function(index_file = file.path("data", "simulated", "simulation_index.csv"),
                                                  model_specs_arg,
                                                  id_var = "id",
                                                  y_var = "y",
                                                  out_dir = file.path("data", "simulated"),
                                                  overwrite = TRUE) {
  idx <- utils::read.csv(index_file, stringsAsFactors = FALSE)
  all_model_vars <- .model_vars_from_specs_cc(model_specs_arg)

  idx$file_cc <- NA_character_
  idx$n_cc <- NA_integer_
  idx$complete_case_rate_common <- NA_real_

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (rr in seq_len(nrow(idx))) {
    mis_file <- .resolve_index_path_cc(idx$file_mis[rr], index_file)
    dat <- utils::read.csv(mis_file, stringsAsFactors = FALSE)

    cc_vars <- unique(c(y_var, id_var, all_model_vars))
    cc_vars <- cc_vars[cc_vars %in% names(dat)]
    cc_rows <- stats::complete.cases(dat[, cc_vars, drop = FALSE])
    dat_cc <- dat[cc_rows, , drop = FALSE]

    base <- tools::file_path_sans_ext(basename(mis_file))
    cc_file <- file.path(out_dir, paste0(base, "_common_cc.csv"))
    if (overwrite || !file.exists(cc_file)) {
      utils::write.csv(dat_cc, cc_file, row.names = FALSE)
    }
    idx$file_cc[rr] <- cc_file
    idx$n_cc[rr] <- nrow(dat_cc)
    idx$complete_case_rate_common[rr] <- nrow(dat_cc) / nrow(dat)

    message("Exported common CC | rep=", idx$rep_id[rr], " | ", idx$mechanism[rr],
            " | n_cc=", nrow(dat_cc), " / ", nrow(dat))
  }

  utils::write.csv(idx, index_file, row.names = FALSE)
  invisible(idx)
}
