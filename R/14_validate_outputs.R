validate_outputs <- function(root) {
  required_files <- c(
    "00_metadata/workflow_dictionary.csv",
    "00_metadata/realized_missingness_rates.csv",
    "00_metadata/complete_case_rates.csv",
    "01_parameter_recovery/estimates_by_replication_all_models.csv",
    "01_parameter_recovery/estimates_by_replication_M2.csv",
    "01_parameter_recovery/param_recovery_all_models.csv",
    "01_parameter_recovery/param_recovery_M2_all_mechanisms.csv",
    "02_model_selection/selected_models_by_replication.csv",
    "02_model_selection/model_selection_rates_all_mechanisms.csv",
    "03_evidence/evidence_by_replication.csv",
    "05_sensitivity/fixed_point_sensitivity.csv",
    "07_runtime/runtime_by_workflow_mechanism.csv",
    "07_runtime/runtime_by_model.csv"
  )
  
  missing <- required_files[!file.exists(file.path(root, required_files))]
  
  if (length(missing) > 0L) {
    stop(
      "The following required output files are missing:\n",
      paste(missing, collapse = "\n")
    )
  }
  
  est_M2 <- readr::read_csv(
    file.path(root, "01_parameter_recovery", "estimates_by_replication_M2.csv"),
    show_col_types = FALSE
  )
  
  if (!all(unique(est_M2$model) == "M2_true")) {
    stop("estimates_by_replication_M2.csv contains models other than M2_true.")
  }
  
  if (any(is.na(est_M2$workflow) | est_M2$workflow == "")) {
    stop("estimates_by_replication_M2.csv contains empty workflow labels.")
  }
  
  message("Output validation passed: ", root)
  invisible(TRUE)
}