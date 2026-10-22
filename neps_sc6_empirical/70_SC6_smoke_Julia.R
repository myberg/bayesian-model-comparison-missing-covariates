# ------------------------------------------------------------
# Julia smoke test from RStudio: robust against spaces in paths
# ------------------------------------------------------------

julia_bin <- Sys.which("julia")
if (julia_bin == "") {
  stop("Julia wurde nicht gefunden. Bitte Julia zum PATH hinzufügen oder julia_bin manuell setzen.")
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Wir arbeiten ab jetzt bewusst mit relativen Pfaden,
# damit der OneDrive-Pfad mit Leerzeichen nicht an Julia übergeben wird.
julia_script_rel <- file.path("neps_sc6_empirical", "71_SC6_run_bayes_medoid.jl")

if (!file.exists(file.path(project_dir, julia_script_rel))) {
  stop("Julia runner not found: ", file.path(project_dir, julia_script_rel))
}

input_files_rel <- list.files(
  path = project_dir,
  pattern = "sc6_wave16_bayes_input\\.csv$",
  recursive = TRUE,
  full.names = FALSE
)

if (length(input_files_rel) == 0L) {
  stop(
    "Keine sc6_wave16_bayes_input.csv gefunden. ",
    "Bitte zuerst RUN_BAYES_ESTIMATION <- TRUE im Master laufen lassen."
  )
}

# Neueste Input-Datei wählen
input_abs <- file.path(project_dir, input_files_rel)
input_info <- file.info(input_abs)
input_rel <- input_files_rel[which.max(input_info$mtime)]

message("Using input file: ", input_rel)

# input_rel sieht ungefähr so aus:
# results/sc6_empirical/sc6_.../04_bayes_cartbb/00_input/sc6_wave16_bayes_input.csv

bayes_input_dir_rel <- dirname(input_rel)
bayes_dir_rel <- dirname(bayes_input_dir_rel)
outdir_rel <- file.path(bayes_dir_rel, "smoke")

dir.create(file.path(project_dir, outdir_rel), recursive = TRUE, showWarnings = FALSE)

message("Using output dir: ", outdir_rel)

# Wichtig: KEIN absoluter Modulpfad und kein --module.
# Das Julia-Skript findet src/BayesREProbit.jl relativ zu sich selbst.
args <- c(
  "--project=.",
  julia_script_rel,
  "--data", input_rel,
  "--outdir", outdir_rel,
  "--medoid_n_iter", "200",
  "--medoid_burnin", "50",
  "--medoid_thin", "10",
  "--n_iter", "300",
  "--burnin", "100",
  "--thin", "5",
  "--chib_reduced_iter", "200",
  "--chib_reduced_burnin", "50",
  "--chib_reduced_thin", "2",
  "--cart_every", "10",
  "--cart_max_depth", "3",
  "--cart_min_samples_leaf", "20",
  "--fixed_points", "mean",
  "--verbose_every", "50"
)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(project_dir)

status <- system2(
  command = julia_bin,
  args = args,
  stdout = "",
  stderr = ""
)

message("Julia exit status: ", status)

list.files(outdir_rel)
