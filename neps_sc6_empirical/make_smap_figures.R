###############################################################################
# make_smap_figures_v4.R
#
# Figures for the SMAP of:
# Bayesian model comparison for random-effects probit models with missing covariates
#
# Produces:
#   Figure 1  Workflow diagram                     -> TikZ standalone PDF
#   Figure 2  Simulation selection performance      -> ggplot from simulation CSV
#   Figure 3  Empirical ranking contrast            -> ggplot from SC6 output CSVs
#   Figure A1 MCMC mixing diagnostics               -> ggplot from draws/diagnostics CSV
#   Figure A2 Runtime/computational burden          -> ggplot from runtime CSV
#
# Usage from project root:
#   source("make_smap_figures.R")
#
# Optional environment variables:
#   FIGURE_ROOT   = project root, default getwd()
#   FIGURE_DIR    = output figure dir, default file.path(FIGURE_ROOT, "figures")
#   RESULTS_DIR   = result root, default file.path(FIGURE_ROOT, "results")
#   SC6_ML_DIR    = optional exact ML benchmark output directory for Figure 3
#   SC6_BAYES_DIR = optional exact final Bayesian output directory for Figure 3/A2
#
# Version 4 changes:
#   - Figure 3 avoids bayes_sensitivity/smoke folders unless SC6_BAYES_DIR is set.
#   - Figure A1 accepts broader diagnostic-column names and no longer stops the full run.
#   - Figure A2 fallback syntax fixed.
#   - Fixed R string escaping in the R-hat diagnostic regex.
#   - The script defines functions only; call make_all_figures() manually.
###############################################################################

required <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr")
missing_required <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0L) {
  stop("Please install required package(s): ", paste(missing_required, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

ROOT <- normalizePath(Sys.getenv("FIGURE_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = file.path(ROOT, "results"))
FIGURE_DIR <- Sys.getenv("FIGURE_DIR", unset = file.path(ROOT, "figures"))
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
TEMPLATE_DIR <- file.path(FIGURE_DIR, "_templates")
dir.create(TEMPLATE_DIR, recursive = TRUE, showWarnings = FALSE)

message_ts <- function(...) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

find_latest_file <- function(pattern,
                             roots = c(RESULTS_DIR, ROOT),
                             exclude = c("smoke", "old", "archive"),
                             prefer_dir = NULL) {
  # Prefer an explicitly supplied directory, e.g. SC6_BAYES_DIR for the final
  # empirical Bayesian run. This avoids accidentally picking a sensitivity run.
  if (!is.null(prefer_dir) && nzchar(prefer_dir) && dir.exists(prefer_dir)) {
    pref <- list.files(prefer_dir, pattern = pattern, recursive = TRUE,
                       full.names = TRUE, ignore.case = TRUE)
    pref <- unique(pref[file.exists(pref)])
    if (length(pref) > 0L) return(pref[which.max(file.info(pref)$mtime)])
  }
  
  roots <- roots[dir.exists(roots)]
  if (length(roots) == 0L) return(NA_character_)
  files <- unlist(lapply(roots, function(r) {
    list.files(r, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  files <- unique(files[file.exists(files)])
  if (length(files) == 0L) return(NA_character_)
  if (length(exclude) > 0L) {
    drop <- Reduce(`|`, lapply(exclude, function(x) grepl(x, files, ignore.case = TRUE)))
    files <- files[!drop]
  }
  if (length(files) == 0L) return(NA_character_)
  files[which.max(file.info(files)$mtime)]
}

read_latest_csv <- function(pattern,
                            roots = c(RESULTS_DIR, ROOT),
                            label = pattern,
                            required = FALSE,
                            exclude = c("smoke", "old", "archive"),
                            prefer_dir = NULL) {
  f <- find_latest_file(pattern, roots = roots, exclude = exclude, prefer_dir = prefer_dir)
  if (is.na(f)) {
    if (required) stop("Could not find required file for ", label, " with pattern: ", pattern, call. = FALSE)
    message_ts("No file found for ", label, ".")
    return(NULL)
  }
  message_ts("Reading ", label, ": ", normalizePath(f, winslash = "/", mustWork = FALSE))
  readr::read_csv(f, show_col_types = FALSE)
}

normalise_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

pick_col_regex <- function(nms, pattern) {
  idx <- grep(pattern, normalise_names(nms), ignore.case = TRUE)
  if (length(idx) == 0L) return(NA_character_)
  nms[idx[1]]
}

save_plot <- function(p, filename, width = 7.2, height = 4.6) {
  pdf_file <- file.path(FIGURE_DIR, paste0(filename, ".pdf"))
  png_file <- file.path(FIGURE_DIR, paste0(filename, ".png"))
  ggplot2::ggsave(pdf_file, p, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(png_file, p, width = width, height = height, dpi = 400)
  message_ts("Wrote ", pdf_file)
  invisible(list(pdf = pdf_file, png = png_file))
}

clean_model_label <- function(x) {
  x <- as.character(x)
  dplyr::recode(x,
                "M1" = "M1: baseline",
                "M2" = "M2: human capital",
                "M3" = "M3: labour history",
                "M4" = "M4: social background",
                "M5" = "M5: full",
                .default = x
  )
}

paper_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_rect(fill = NA),
      plot.title.position = "plot"
    )
}

# -----------------------------------------------------------------------------
# Figure 1: Workflow diagram via TikZ standalone
# -----------------------------------------------------------------------------

write_workflow_tikz <- function(out_tex = file.path(FIGURE_DIR, "fig01_workflow_tikz_standalone.tex")) {
  tikz <- '\\documentclass[tikz,border=6pt]{standalone}
\\usepackage{amsmath}
\\usepackage{tikz}
\\usetikzlibrary{arrows.meta,positioning,fit,calc,shapes.geometric}

\\begin{document}
\\begin{tikzpicture}[
  font=\\sffamily\\small,
  node distance=8mm and 10mm,
  >=Latex,
  box/.style={draw, rounded corners=2pt, align=center, inner sep=5pt, minimum height=9mm, text width=34mm},
  widebox/.style={draw, rounded corners=2pt, align=center, inner sep=5pt, minimum height=9mm, text width=45mm},
  target/.style={draw, rounded corners=2pt, align=center, inner sep=5pt, minimum height=9mm, text width=42mm, thick},
  group/.style={draw, rounded corners=3pt, inner sep=6pt},
  arrow/.style={->, thick}
]

\\node[widebox] (data) {Observed survey data\\\\[-1mm]
$y$ observed, $X=(X_{\\mathrm{obs}},X_{\\mathrm{mis}})$};

\\node[box, below left=15mm and 26mm of data] (mice) {Common MICE--CART\\\\imputation model};
\\node[box, below=of mice] (glmm) {Random-intercept\\\\probit GLMMs\\\\$M_1,\\ldots,M_5$};
\\node[target, below=of glmm] (mlcrit) {ML benchmark\\\\AIC, BIC, pooled LRT,\\\\BIC-based log BF};

\\node[box, below right=15mm and 26mm of data] (cartbb) {CART--BB reference chain\\\\on union of covariates};
\\node[box, below=of cartbb] (medoid) {Select one common\\\\medoid completion\\\\$X_{\\mathrm{med}}$};
\\node[box, below=of medoid] (gibbs) {Fixed-design Gibbs runs\\\\$M_1,\\ldots,M_5$ on same $X_{\\mathrm{med}}$};
\\node[target, below=of gibbs] (bayescrit) {Bayesian workflow\\\\conditional Chib evidence,\\\\WAIC, sensitivity checks};

\\draw[arrow] (data) -- (mice);
\\draw[arrow] (mice) -- (glmm);
\\draw[arrow] (glmm) -- (mlcrit);
\\draw[arrow] (data) -- (cartbb);
\\draw[arrow] (cartbb) -- (medoid);
\\draw[arrow] (medoid) -- (gibbs);
\\draw[arrow] (gibbs) -- (bayescrit);

\\node[widebox, below=16mm of $(mlcrit)!0.5!(bayescrit)$] (compare) {Compare model rankings and evidence differences\\\\not causal effects of individual covariates};
\\draw[arrow] (mlcrit) |- (compare);
\\draw[arrow] (bayescrit) |- (compare);

\\node[group, fit=(mice)(glmm)(mlcrit), label={[font=\\sffamily\\bfseries]above:ML/MICE benchmark}] {};
\\node[group, fit=(cartbb)(medoid)(gibbs)(bayescrit), label={[font=\\sffamily\\bfseries]above:Bayesian conditional-evidence workflow}] {};

\\node[align=center, font=\\sffamily\\footnotesize, below=2mm of bayescrit, text width=45mm] {Target: $p(y\\mid X_{\\mathrm{med}},M_k)$,\\\\not observed-data BF integrating over $X_{\\mathrm{mis}}$.};

\\end{tikzpicture}
\\end{document}
'
writeLines(tikz, out_tex)
message_ts("Wrote ", out_tex)
invisible(out_tex)
}

compile_workflow_tikz <- function() {
  tex <- write_workflow_tikz()
  pdflatex <- Sys.which("pdflatex")
  if (pdflatex == "") {
    message_ts("pdflatex not found. Compile manually: pdflatex ", basename(tex))
    return(invisible(tex))
  }
  old <- getwd(); on.exit(setwd(old), add = TRUE)
  setwd(FIGURE_DIR)
  status <- system2(pdflatex, args = c("-interaction=nonstopmode", basename(tex)), stdout = TRUE, stderr = TRUE)
  out_pdf <- sub("\\.tex$", ".pdf", tex)
  final_pdf <- file.path(FIGURE_DIR, "fig01_workflow.pdf")
  if (file.exists(out_pdf)) file.copy(out_pdf, final_pdf, overwrite = TRUE)
  message_ts("Wrote ", final_pdf)
  invisible(final_pdf)
}

# -----------------------------------------------------------------------------
# Figure 2: Simulation model-selection performance
# -----------------------------------------------------------------------------

standardize_sim_selection <- function(dat) {
  original_names <- names(dat)
  names(dat) <- stringr::str_replace_all(names(dat), "\\s+", "_")
  norm <- normalise_names(names(dat))
  
  # Case 1: already in long format.
  method_col <- intersect(names(dat), c("method", "Method", "workflow", "Workflow", "ID", "id"))[1]
  crit_col <- intersect(names(dat), c("criterion", "Criterion", "criteria", "Criteria"))[1]
  mech_col <- intersect(names(dat), c("mechanism", "Mechanism", "missingness", "Missingness", "scenario", "Scenario"))[1]
  rate_col <- intersect(names(dat), c("selection_rate", "Selection_rate", "selection_frequency", "Selection_frequency", "selected_rate", "M2_selected"))[1]
  mcse_col <- intersect(names(dat), c("mcse", "MCSE", "se", "SE"))[1]
  
  if (!is.na(method_col) && !is.na(crit_col) && !is.na(mech_col) && !is.na(rate_col)) {
    out <- dat %>%
      rename(method = all_of(method_col), criterion = all_of(crit_col),
             mechanism = all_of(mech_col), selection_rate = all_of(rate_col)) %>%
      mutate(
        method = as.character(method),
        criterion = as.character(criterion),
        mechanism = as.character(mechanism),
        selection_rate = as.numeric(selection_rate),
        mcse = if (!is.na(mcse_col)) as.numeric(.data[[mcse_col]]) else NA_real_
      ) %>%
      select(method, criterion, mechanism, selection_rate, mcse)
    return(out)
  }
  
  # Case 2: paper-table wide format with ID/criterion and MCAR/MAR/MNAR columns.
  # We use only columns that look like selection frequencies, not log-BF columns.
  id_col <- intersect(names(dat), c("ID", "id", "method", "Method", "workflow", "Workflow"))[1]
  crit_col <- intersect(names(dat), c("Criterion", "criterion", "criteria", "Criteria"))[1]
  if (is.na(id_col)) id_col <- names(dat)[1]
  if (is.na(crit_col) && ncol(dat) >= 2) crit_col <- names(dat)[2]
  
  candidate_mech <- names(dat)[grepl("mcar|mar|mnar", normalise_names(names(dat)))]
  candidate_mech <- candidate_mech[!grepl("log|bf|mean|sd|se|evidence", normalise_names(candidate_mech))]
  
  if (length(candidate_mech) >= 2L) {
    tmp <- dat %>%
      select(all_of(c(id_col, crit_col, candidate_mech))) %>%
      rename(method = all_of(id_col), criterion = all_of(crit_col)) %>%
      pivot_longer(cols = all_of(candidate_mech), names_to = "mechanism", values_to = "selection_rate") %>%
      mutate(
        method = as.character(method),
        criterion = as.character(criterion),
        mechanism = case_when(
          grepl("mcar", mechanism, ignore.case = TRUE) ~ "MCAR",
          grepl("mnar", mechanism, ignore.case = TRUE) ~ "MNAR",
          grepl("mar", mechanism, ignore.case = TRUE) ~ "MAR-type",
          TRUE ~ mechanism
        ),
        selection_rate = suppressWarnings(as.numeric(selection_rate)),
        mcse = NA_real_
      ) %>%
      filter(is.finite(selection_rate), selection_rate >= 0, selection_rate <= 1)
    if (nrow(tmp) > 0L) return(tmp)
  }
  
  attr(dat, "original_names") <- original_names
  return(NULL)
}

make_fig02_sim_selection <- function() {
  sim_prefer <- Sys.getenv("SIM_SELECTION_FILE", unset = "")
  dat <- NULL
  if (nzchar(sim_prefer) && file.exists(sim_prefer)) {
    message_ts("Reading simulation selection/evidence summary: ", normalizePath(sim_prefer, winslash = "/", mustWork = FALSE))
    dat <- readr::read_csv(sim_prefer, show_col_types = FALSE)
  } else {
    dat <- read_latest_csv(
      pattern = "(sim.*selection.*evidence|simulation.*selection|selection.*evidence|selection_accuracy|model_selection).*\\.csv$",
      roots = c(RESULTS_DIR, file.path(ROOT, "outputs"), ROOT),
      exclude = c("smoke", "old", "archive", "_templates"),
      label = "simulation selection/evidence summary"
    )
  }
  
  if (is.null(dat)) {
    template <- tibble::tribble(
      ~method, ~criterion, ~mechanism, ~selection_rate, ~mcse,
      "BD-ML", "BIC", "MCAR", NA_real_, NA_real_,
      "CC-ML", "BIC", "MCAR", NA_real_, NA_real_,
      "MICE-CART+ML", "BIC", "MCAR", NA_real_, NA_real_,
      "missForest+ML", "BIC", "MCAR", NA_real_, NA_real_,
      "CART-BB+Bayes", "Chib evidence", "MCAR", NA_real_, NA_real_,
      "CART-BB+Bayes", "WAIC", "MCAR", NA_real_, NA_real_,
      "BD-ML", "BIC", "MAR-type", NA_real_, NA_real_,
      "CC-ML", "BIC", "MAR-type", NA_real_, NA_real_,
      "MICE-CART+ML", "BIC", "MAR-type", NA_real_, NA_real_,
      "missForest+ML", "BIC", "MAR-type", NA_real_, NA_real_,
      "CART-BB+Bayes", "Chib evidence", "MAR-type", NA_real_, NA_real_,
      "CART-BB+Bayes", "WAIC", "MAR-type", NA_real_, NA_real_,
      "BD-ML", "BIC", "MNAR", NA_real_, NA_real_,
      "CC-ML", "BIC", "MNAR", NA_real_, NA_real_,
      "MICE-CART+ML", "BIC", "MNAR", NA_real_, NA_real_,
      "missForest+ML", "BIC", "MNAR", NA_real_, NA_real_,
      "CART-BB+Bayes", "Chib evidence", "MNAR", NA_real_, NA_real_,
      "CART-BB+Bayes", "WAIC", "MNAR", NA_real_, NA_real_
    )
    out <- file.path(TEMPLATE_DIR, "fig02_sim_selection_template.csv")
    readr::write_csv(template, out)
    message_ts("No simulation summary found. Fill template or set SIM_SELECTION_FILE and rerun: ", out)
    return(invisible(NULL))
  }
  
  plotdat <- standardize_sim_selection(dat)
  if (is.null(plotdat) || nrow(plotdat) == 0L) {
    out <- file.path(TEMPLATE_DIR, "fig02_sim_selection_unreadable_columns.txt")
    writeLines(c("Could not detect columns for Figure 2.", "Available columns:", paste(names(dat), collapse = ", ")), out)
    message_ts("Simulation file found but columns could not be detected. See: ", out)
    return(invisible(NULL))
  }
  
  plotdat <- plotdat %>%
    mutate(
      criterion = as.character(criterion),
      label = paste(method, criterion, sep = " | "),
      mechanism = factor(mechanism, levels = c("MCAR", "MAR-type", "MAR", "MNAR")),
      lo = ifelse(is.na(mcse), NA_real_, pmax(0, selection_rate - 1.96 * mcse)),
      hi = ifelse(is.na(mcse), NA_real_, pmin(1, selection_rate + 1.96 * mcse))
    ) %>%
    filter(is.finite(selection_rate), !is.na(mechanism))
  
  if (nrow(plotdat) == 0L) {
    message_ts("Simulation selection data contain no finite selection rates.")
    return(invisible(NULL))
  }
  
  p <- ggplot(plotdat, aes(x = selection_rate, y = reorder(label, selection_rate))) +
    geom_point(size = 2.1) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, na.rm = TRUE) +
    facet_wrap(~ mechanism, nrow = 1) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    labs(x = "Selection frequency of the data-generating model", y = NULL) +
    paper_theme(base_size = 9)
  
  save_plot(p, "fig02_sim_selection_performance", width = 8.2, height = 4.8)
}

# -----------------------------------------------------------------------------
# Figure 3: Empirical ranking contrast across criteria
# -----------------------------------------------------------------------------

make_fig03_empirical_ranking <- function() {
  ml_prefer <- Sys.getenv("SC6_ML_DIR", unset = "")
  bayes_prefer <- Sys.getenv("SC6_BAYES_DIR", unset = "")
  
  ml <- read_latest_csv("ml_fit_statistics_summary\\.csv$",
                        label = "ML fit statistics",
                        required = TRUE,
                        prefer_dir = ml_prefer)
  
  # For the main empirical figure, do not accidentally use bayes_sensitivity.
  # Set SC6_BAYES_DIR to the final/paper Bayesian output directory if needed.
  bayes_results <- read_latest_csv("bayes_model_results\\.csv$",
                                   label = "Bayesian model results",
                                   exclude = c("smoke", "old", "archive", "bayes_sensitivity"),
                                   prefer_dir = bayes_prefer)
  bayes_rank <- read_latest_csv("bayes_conditional_logml_ranking\\.csv$",
                                label = "Bayesian conditional logML ranking",
                                exclude = c("smoke", "old", "archive", "bayes_sensitivity"),
                                prefer_dir = bayes_prefer)
  
  ml_plot <- ml %>%
    transmute(
      model = model,
      AIC = if ("AIC_median" %in% names(ml)) AIC_median else AIC,
      BIC = if ("BIC_median" %in% names(ml)) BIC_median else BIC
    ) %>%
    pivot_longer(c(AIC, BIC), names_to = "criterion", values_to = "value") %>%
    group_by(criterion) %>%
    mutate(loss = value - min(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(criterion = paste0(criterion, " (ML)"))
  
  if (!is.null(bayes_rank) && "log_marginal_likelihood" %in% names(bayes_rank)) {
    br <- bayes_rank
    if ("fixed_point" %in% names(br)) br <- br %>% filter(fixed_point == "mean")
    chib <- br %>%
      transmute(model, criterion = "Chib evidence (Bayes)", value = log_marginal_likelihood) %>%
      group_by(criterion) %>%
      mutate(loss = max(value, na.rm = TRUE) - value) %>%
      ungroup()
  } else if (!is.null(bayes_results) && "log_marginal_likelihood" %in% names(bayes_results)) {
    chib <- bayes_results %>%
      transmute(model, criterion = "Chib evidence (Bayes)", value = log_marginal_likelihood) %>%
      group_by(criterion) %>%
      mutate(loss = max(value, na.rm = TRUE) - value) %>%
      ungroup()
  } else {
    stop("Could not find Bayesian log marginal likelihoods for Figure 3. Set SC6_BAYES_DIR to the final Bayesian output directory.", call. = FALSE)
  }
  
  if (!is.null(bayes_results) && "WAIC" %in% names(bayes_results)) {
    waic <- bayes_results %>%
      transmute(model, criterion = "WAIC (Bayes)", value = WAIC) %>%
      group_by(criterion) %>%
      mutate(loss = value - min(value, na.rm = TRUE)) %>%
      ungroup()
  } else {
    waic <- NULL
  }
  
  plotdat <- bind_rows(ml_plot, chib, waic) %>%
    mutate(
      criterion = factor(criterion, levels = c("AIC (ML)", "BIC (ML)", "Chib evidence (Bayes)", "WAIC (Bayes)")),
      model_label = factor(clean_model_label(model), levels = rev(clean_model_label(paste0("M", 1:5)))),
      best = loss == 0
    )
  
  p <- ggplot(plotdat, aes(x = loss, y = model_label)) +
    geom_segment(aes(x = 0, xend = loss, yend = model_label), linewidth = 0.35, na.rm = TRUE) +
    geom_point(aes(shape = best), size = 2.4, na.rm = TRUE) +
    facet_wrap(~ criterion, scales = "free_x", nrow = 1) +
    labs(x = "Loss relative to best model (criterion-specific scale)", y = NULL, shape = "Best") +
    paper_theme()
  
  save_plot(p, "fig03_empirical_ranking_contrast", width = 8.0, height = 4.0)
}

# -----------------------------------------------------------------------------
# Figure A1: MCMC mixing diagnostics
# -----------------------------------------------------------------------------

make_figA1_mcmc_diagnostics <- function() {
  diag <- read_latest_csv("(bayes_)?mcmc_diagnostics.*\\.csv$|bayes_diagnostics.*\\.csv$|table_mcmc_diagnostics_compact\\.csv$",
                          label = "Bayesian MCMC diagnostics")
  draws <- read_latest_csv("(bayes_)?draws_long.*\\.csv$|posterior_draws_long.*\\.csv$",
                           label = "Bayesian posterior draws long")
  
  if (!is.null(draws)) {
    names(draws) <- stringr::str_replace_all(names(draws), "\\s+", "_")
    model_col <- intersect(names(draws), c("model", "Model"))[1]
    draw_col <- intersect(names(draws), c("draw", "iter", "iteration", "Iteration"))[1]
    par_col <- intersect(names(draws), c("parameter", "term", "variable", "name"))[1]
    val_col <- intersect(names(draws), c("value", "draw_value", "estimate"))[1]
    if (any(is.na(c(model_col, draw_col, par_col, val_col)))) {
      message_ts("Posterior-draws file found but required columns could not be detected. Expected model, draw, parameter, value. Trying diagnostics file instead.")
    } else {
      key_terms <- c("sigma2_a", "logLik", "log_likelihood", "age_z", "female", "years_education", "unemp_ge5")
      plotdat <- draws %>%
        rename(model = all_of(model_col), draw = all_of(draw_col), parameter = all_of(par_col), value = all_of(val_col)) %>%
        mutate(parameter = as.character(parameter)) %>%
        filter(model %in% c("M3", "M5"), parameter %in% key_terms) %>%
        mutate(panel = paste(model, parameter, sep = ": "))
      if (nrow(plotdat) == 0L) {
        plotdat <- draws %>%
          rename(model = all_of(model_col), draw = all_of(draw_col), parameter = all_of(par_col), value = all_of(val_col)) %>%
          filter(model %in% c("M3", "M5")) %>%
          group_by(model, parameter) %>%
          mutate(n_par = cur_group_id()) %>%
          ungroup() %>%
          filter(n_par <= 8) %>%
          mutate(panel = paste(model, parameter, sep = ": "))
      }
      p <- ggplot(plotdat, aes(x = as.numeric(draw), y = as.numeric(value))) +
        geom_line(linewidth = 0.2, na.rm = TRUE) +
        facet_wrap(~ panel, scales = "free_y") +
        labs(x = "Retained draw", y = "Parameter value") +
        paper_theme(base_size = 9)
      save_plot(p, "figA1_mcmc_mixing_diagnostics", width = 7.4, height = 5.2)
      return(invisible(TRUE))
    }
  }
  
  if (!is.null(diag)) {
    original_names <- names(diag)
    names(diag) <- stringr::str_replace_all(names(diag), "\\s+", "_")
    norm <- normalise_names(names(diag))
    
    model_col <- pick_col_regex(names(diag), "^(model|m)$|candidate_model")
    if (is.na(model_col)) model_col <- pick_col_regex(names(diag), "method|workflow")
    par_col <- pick_col_regex(names(diag), "parameter|term|variable|name|estimand|quantity")
    
    ess_cols <- names(diag)[grep("ess|effective_sample", norm)]
    rhat_cols <- names(diag)[grep("rhat|r_hat|r\\.hat", norm)]
    mcse_cols <- names(diag)[grep("mcse|monte_carlo|se_mc", norm)]
    keep_cols <- unique(c(ess_cols, rhat_cols, mcse_cols))
    
    if (length(keep_cols) == 0L) {
      out <- file.path(TEMPLATE_DIR, "figA1_mcmc_diagnostics_columns_found.txt")
      writeLines(c("Diagnostics file found, but no ESS/Rhat/MCSE columns were detected.",
                   "Columns found:", paste(original_names, collapse = ", ")), out)
      message_ts("Diagnostics file found but no ESS/Rhat/MCSE columns were detected. Column list written to: ", out)
      return(invisible(NULL))
    }
    
    tmp <- diag
    if (is.na(model_col)) tmp$model <- "All models/workflows" else tmp$model <- as.character(tmp[[model_col]])
    if (is.na(par_col)) {
      tmp$parameter <- seq_len(nrow(tmp))
      tmp$parameter <- paste0("row ", tmp$parameter)
    } else {
      tmp$parameter <- as.character(tmp[[par_col]])
    }
    
    plotdat <- tmp %>%
      select(model, parameter, all_of(keep_cols)) %>%
      pivot_longer(cols = all_of(keep_cols), names_to = "diagnostic", values_to = "value") %>%
      mutate(
        value = suppressWarnings(as.numeric(value)),
        diagnostic_clean = case_when(
          grepl("ess|effective", diagnostic, ignore.case = TRUE) ~ "Effective sample size",
          grepl("rhat|r_hat|r\\.hat", diagnostic, ignore.case = TRUE) ~ "R-hat",
          grepl("mcse|monte", diagnostic, ignore.case = TRUE) ~ "MCSE",
          TRUE ~ diagnostic
        )
      ) %>%
      filter(is.finite(value))
    
    # Prefer empirical M3/M5 if present; otherwise use all rows in the compact diagnostic table.
    if (any(plotdat$model %in% c("M3", "M5"))) {
      plotdat <- plotdat %>% filter(model %in% c("M3", "M5"))
    }
    
    if (nrow(plotdat) == 0L) {
      message_ts("Diagnostics file was parsed, but no finite diagnostic values were available for plotting.")
      return(invisible(NULL))
    }
    
    # Keep the plot readable if the compact table has many rows/diagnostics.
    plotdat <- plotdat %>%
      group_by(diagnostic_clean, model) %>%
      slice_head(n = 25) %>%
      ungroup()
    
    p <- ggplot(plotdat, aes(x = value, y = reorder(parameter, value))) +
      geom_point(size = 1.8, na.rm = TRUE) +
      facet_grid(diagnostic_clean ~ model, scales = "free", space = "free") +
      labs(x = NULL, y = NULL) +
      paper_theme(base_size = 9)
    save_plot(p, "figA1_mcmc_mixing_diagnostics", width = 7.4, height = 5.2)
    return(invisible(TRUE))
  }
  
  template <- tibble::tribble(
    ~model, ~parameter, ~ess_bulk, ~rhat, ~mcse,
    "M3", "sigma2_a", NA_real_, NA_real_, NA_real_,
    "M3", "unemp_ge5", NA_real_, NA_real_, NA_real_,
    "M5", "sigma2_a", NA_real_, NA_real_, NA_real_,
    "M5", "unemp_ge5", NA_real_, NA_real_, NA_real_
  )
  out <- file.path(TEMPLATE_DIR, "figA1_mcmc_diagnostics_template.csv")
  readr::write_csv(template, out)
  message_ts("No MCMC diagnostics or posterior-draws file found. Fill template and rerun: ", out)
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Figure A2: Runtime and computational burden
# -----------------------------------------------------------------------------

standardize_runtime <- function(runtime) {
  original_names <- names(runtime)
  names(runtime) <- stringr::str_replace_all(names(runtime), "\\s+", "_")
  nms <- names(runtime)
  norm <- normalise_names(nms)
  
  pick <- function(pattern) {
    idx <- grep(pattern, norm, ignore.case = TRUE)
    if (length(idx) == 0L) return(NA_character_)
    nms[idx[1]]
  }
  
  workflow_col <- pick("^(workflow|method|approach|algorithm|estimator|procedure|stage)$|workflow|method|approach|algorithm|estimator|procedure|stage")
  mech_col <- pick("mechanism|missingness|scenario|condition")
  model_col <- pick("^model$|model")
  time_col <- pick("runtime_seconds|seconds|elapsed_seconds|time_seconds|runtime_sec|run_time_sec|elapsed_sec|duration_sec|wall.*sec|time_sec")
  
  minutes_col <- pick("runtime_minutes|minutes|elapsed_minutes|time_minutes|duration_min|wall.*min")
  hours_col <- pick("runtime_hours|hours|elapsed_hours|time_hours|duration_h")
  
  if (is.na(time_col) && is.na(minutes_col) && is.na(hours_col)) {
    # Last-resort: use a numeric column with runtime/time/elapsed in its name.
    cand <- nms[grepl("runtime|run_time|elapsed|duration|time", norm)]
    cand <- cand[vapply(runtime[cand], function(x) is.numeric(x) || is.integer(x), logical(1))]
    if (length(cand) > 0L) time_col <- cand[1]
  }
  
  if (is.na(workflow_col)) {
    # Last-resort: create a workflow column from the first character-like column.
    char_cols <- nms[vapply(runtime, function(x) is.character(x) || is.factor(x), logical(1))]
    workflow_col <- if (length(char_cols) > 0L) char_cols[1] else NA_character_
  }
  
  if (is.na(workflow_col) || (is.na(time_col) && is.na(minutes_col) && is.na(hours_col))) {
    attr(runtime, "original_names") <- original_names
    return(NULL)
  }
  
  out <- runtime %>%
    mutate(
      workflow = as.character(.data[[workflow_col]]),
      mechanism = if (!is.na(mech_col)) as.character(.data[[mech_col]]) else "All/empirical",
      model = if (!is.na(model_col)) as.character(.data[[model_col]]) else "all",
      runtime_minutes = case_when(
        !is.na(time_col) ~ suppressWarnings(as.numeric(.data[[time_col]]) / 60),
        !is.na(minutes_col) ~ suppressWarnings(as.numeric(.data[[minutes_col]])),
        !is.na(hours_col) ~ suppressWarnings(as.numeric(.data[[hours_col]]) * 60),
        TRUE ~ NA_real_
      )
    ) %>%
    select(workflow, mechanism, model, runtime_minutes) %>%
    filter(is.finite(runtime_minutes), runtime_minutes > 0)
  
  out
}

make_figA2_runtime <- function() {
  runtime_prefer <- Sys.getenv("RUNTIME_FILE", unset = "")
  runtime <- NULL
  if (nzchar(runtime_prefer) && file.exists(runtime_prefer)) {
    message_ts("Reading runtime summary: ", normalizePath(runtime_prefer, winslash = "/", mustWork = FALSE))
    runtime <- readr::read_csv(runtime_prefer, show_col_types = FALSE)
  } else {
    runtime <- read_latest_csv("(runtime_summary|simulation_runtime|benchmark_runtime|run_time).*\\.csv$",
                               roots = c(RESULTS_DIR, file.path(ROOT, "outputs"), ROOT),
                               exclude = c("smoke", "old", "archive", "_templates"),
                               label = "runtime summary")
  }
  
  plotdat <- NULL
  if (!is.null(runtime)) {
    plotdat <- standardize_runtime(runtime)
    if (is.null(plotdat) || nrow(plotdat) == 0L) {
      out <- file.path(TEMPLATE_DIR, "figA2_runtime_unreadable_columns.txt")
      writeLines(c("Could not detect columns for Figure A2.", "Available columns:", paste(names(runtime), collapse = ", ")), out)
      message_ts("Runtime file found but required columns could not be detected. See: ", out)
      runtime <- NULL
    }
  }
  
  if (is.null(runtime)) {
    # Empirical fallback: model-specific Bayesian runtimes, not a full simulation benchmark.
    bayes_prefer <- Sys.getenv("SC6_BAYES_DIR", unset = "")
    bayes_results <- read_latest_csv("bayes_model_results\\.csv$",
                                     label = "Bayesian model results for runtime fallback",
                                     exclude = c("smoke", "old", "archive", "bayes_sensitivity"),
                                     prefer_dir = bayes_prefer)
    if (!is.null(bayes_results) && "runtime_seconds" %in% names(bayes_results)) {
      plotdat <- bayes_results %>%
        transmute(workflow = "Bayesian fixed-medoid Gibbs",
                  mechanism = "NEPS SC6",
                  model = model,
                  runtime_minutes = as.numeric(runtime_seconds) / 60) %>%
        filter(is.finite(runtime_minutes), runtime_minutes > 0)
    }
  }
  
  if (is.null(plotdat) || nrow(plotdat) == 0L) {
    template <- tibble::tribble(
      ~workflow, ~mechanism, ~model, ~runtime_seconds,
      "BD-ML", "MCAR", "all", NA_real_,
      "CC-ML", "MCAR", "all", NA_real_,
      "MICE-CART+ML", "MCAR", "all", NA_real_,
      "missForest+ML", "MCAR", "all", NA_real_,
      "CART-BB reference chain", "MCAR", "reference", NA_real_,
      "Bayesian fixed-medoid Gibbs", "MCAR", "all", NA_real_,
      "Reduced-Gibbs Chib", "MCAR", "all", NA_real_,
      "BD-ML", "MAR-type", "all", NA_real_,
      "CC-ML", "MAR-type", "all", NA_real_,
      "MICE-CART+ML", "MAR-type", "all", NA_real_,
      "missForest+ML", "MAR-type", "all", NA_real_,
      "CART-BB reference chain", "MAR-type", "reference", NA_real_,
      "Bayesian fixed-medoid Gibbs", "MAR-type", "all", NA_real_,
      "Reduced-Gibbs Chib", "MAR-type", "all", NA_real_,
      "BD-ML", "MNAR", "all", NA_real_,
      "CC-ML", "MNAR", "all", NA_real_,
      "MICE-CART+ML", "MNAR", "all", NA_real_,
      "missForest+ML", "MNAR", "all", NA_real_,
      "CART-BB reference chain", "MNAR", "reference", NA_real_,
      "Bayesian fixed-medoid Gibbs", "MNAR", "all", NA_real_,
      "Reduced-Gibbs Chib", "MNAR", "all", NA_real_
    )
    out <- file.path(TEMPLATE_DIR, "figA2_runtime_template.csv")
    readr::write_csv(template, out)
    message_ts("No usable runtime data found. Fill template or set RUNTIME_FILE and rerun: ", out)
    return(invisible(NULL))
  }
  
  plotdat <- plotdat %>%
    group_by(workflow, mechanism) %>%
    summarise(runtime_minutes = median(runtime_minutes, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(plotdat, aes(x = runtime_minutes, y = reorder(workflow, runtime_minutes))) +
    geom_point(size = 2.0) +
    facet_wrap(~ mechanism, scales = "free_x") +
    scale_x_log10() +
    labs(x = "Runtime in minutes (log scale)", y = NULL) +
    paper_theme(base_size = 9)
  
  save_plot(p, "figA2_runtime_computational_burden", width = 7.4, height = 4.8)
}

# Quick source inspection helper. This is useful when Figure 2 or A2 is missing.
inspect_figure_sources <- function(max_files = 30L) {
  roots <- unique(c(RESULTS_DIR, file.path(ROOT, "outputs"), ROOT))
  roots <- roots[dir.exists(roots)]
  csvs <- unique(unlist(lapply(roots, function(r) list.files(r, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)), use.names = FALSE))
  csvs <- csvs[!grepl("_templates|smoke|archive", csvs, ignore.case = TRUE)]
  cand <- csvs[grepl("selection|evidence|runtime|run_time|elapsed|diagnostic", basename(csvs), ignore.case = TRUE)]
  cand <- head(cand[order(file.info(cand)$mtime, decreasing = TRUE)], max_files)
  if (length(cand) == 0L) {
    message_ts("No candidate figure-source CSVs found.")
    return(invisible(NULL))
  }
  out <- lapply(cand, function(f) {
    d <- tryCatch(readr::read_csv(f, n_max = 1, show_col_types = FALSE), error = function(e) NULL)
    tibble::tibble(file = normalizePath(f, winslash = "/", mustWork = FALSE),
                   columns = if (is.null(d)) "<read error>" else paste(names(d), collapse = ", "))
  }) %>% bind_rows()
  print(out, n = nrow(out), width = Inf)
  invisible(out)
}

# -----------------------------------------------------------------------------
# Run all figures
# -----------------------------------------------------------------------------

make_all_figures <- function() {
  compile_workflow_tikz()
  make_fig02_sim_selection()
  make_fig03_empirical_ranking()
  make_figA1_mcmc_diagnostics()
  make_figA2_runtime()
  message_ts("Done. Figures are in: ", normalizePath(FIGURE_DIR, winslash = "/", mustWork = FALSE))
}

# Deliberately do not run automatically. Use:
#   source("make_smap_figures_v2.R")
#   make_all_figures()
