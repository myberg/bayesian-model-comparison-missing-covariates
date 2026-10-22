# Execution guide and retained script map

All paths below are relative to the repository root. This package retains the source-level computational pathway; it has **not** passed a complete R/Julia execution test in the preparation environment.

## 1. Simulation analysis

Main driver: `R/execute_master_v14.R`.

Its first section loads `00_config.R`, the data-generating process, CART-BB and Bayesian samplers, ML/MICE benchmark functions, the Julia bridge, necessary model-comparison helpers and the table-generating functions. A subsequent section loads the result harmonization, metadata, diagnostics and paper tables. Two optional standalone drivers are `R/execute_ml_only.R` and `R/run_prior_sensitivity_MED_Bayes.R`.

`RUN_MODE` is **set to `full100` inside the master**. Output directories may be recreated. Change this to `smoke1` for first local trials **in a disposable workspace**, and check all directories and options before longer runs. The source files have not been normalized to a new execution standard because that could alter the computational history.

## 2. Julia simulation implementation

The shared numerical implementation is `julia/src/BayesREProbit.jl`; `julia/run_bayes.jl` supplies the entry point. `julia/src/diagnostic_draw_export.jl` defines optional diagnostic exports. The `julia/Project.toml` and `julia/Manifest.toml` come from the later working environment, with the package name normalized. Julia may need to refresh its manifest metadata when resolving the project.

For a first local check:

```powershell
julia --project=julia -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project=julia julia/test/conditional_medoid_smoke.jl
```

The full package/runtime check is pending.

## 3. R tests

```powershell
Rscript R/test_conditional_chib_smoke.R
Rscript R/test_ml_d2_smoke.R
Rscript R/test_ml_d3_smoke.R
```

The ML tests now reference the existing `R/09_ml_imputation_estimates.R`.
Their successful execution still has to be confirmed in a suitable R environment.

## 4. NEPS SC6 empirical analysis

These scripts must run with access-controlled NEPS inputs and output paths configured as in `.env.example`:

1. `50_Datenaufbereitung_SC6_paper_ready.R` – derive analysis dataset
2. `51_SC6_data_checks.R` and `52_SC6_descriptives.R` – validation and descriptive statistics
3. `60_SC6_ml_mice_cart_glmer.R` – ML and MICE comparisons
4. `70_SC6_bayes_cartbb_models.R` – Julia input preparation and command construction
5. `71_SC6_run_bayes_medoid.jl` – common medoid and conditional evidence; now loads the shared Julia module at `julia/src/BayesREProbit.jl`
6. `73_SC6_postprocess_bayes_results.R` and `80_SC6_make_paper_tables.R` – output summaries
7. `91_SC6_bayes_sensitivity_driver.R` – extra sensitivity runs
8. `make_smap_figures.R` – paper-oriented figures from saved output data

`70_SC6_smoke_Julia.R` is a small technical trial script that is not a substitute for a full empirical reproduction.
The Julia core previously had an almost identical copy under the empirical directory.

Do not infer historic command-line flags from currently generated shell or batch scripts.
Recorded Bayesian evidence values from the saved June 2026 working run align with key published estimates, but the complete historical execution command is not preserved.
