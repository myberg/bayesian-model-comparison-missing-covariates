# Bayesian model comparison with incomplete covariates

**Research code for:** Bergrab, M. (2026). *Bayesian model comparison for random effects probit models with missing covariates*. **Statistical Methods & Applications**. [https://doi.org/10.1007/s10260-026-00878-5](https://doi.org/10.1007/s10260-026-00878-5)

This repository represents the state of the final published results as a a **curated code archive derived from publication working directory**.
It brings together the R and Julia implementations associated with the published study, while excluding intermediate working files and restricted data.
The published article is the reference for reported findings.

## What is here

| Directory | Contents |
| --- | --- |
| `R/` | Simulation design, missing-data methods, maximum-likelihood and imputation benchmarks, Julia interface, diagnostics, sensitivity analyses, and tables and figures |
| `julia/` | Bayesian random-effects probit implementation, dependency specifications, and a synthetic smoke test |
| `neps_sc6_empirical/` | NEPS Starting Cohort 6, Wave 16: data preparation, empirical analyses, model comparison, and output generation; the Julia runner uses the shared implementation in `julia/src/` |
| `docs/` | Execution instructions, data-access requirements, script inventory, and quality-review notes |
| `tests/` | Static integrity checks that do not require NEPS data, R, or Julia |

The main simulation driver is `R/execute_master.R`. 
The empirical analysis comprises the numbered scripts in `neps_sc6_empirical/`.
See [execution instructions](docs/EXECUTION.md) for the script order and computational requirements.

No NEPS individual records, prepared respondent-level files, simulation exports, execution logs or unpublished results have been placed in this repository.
Historic working data, failed experiments and previous source generations remain only in the private working directory.
Using the NEPS data requires a data protection contract with the Leibniz-Institute for educational trajectories.
For details see [www.neps-data.de](www.neps-data.de).


## Method and interpretation

The workflow addresses random-effects probit models with incompletely observed covariates.
Missing covariates are updated using a decision-tree-assisted Bayesian-bootstrap procedure (CART-BB).
R provides the simulation and comparison workflows, while Julia implements the Bayesian estimation and evidence calculations.
The preserved empirical Chib calculations are conditional on one fixed, common CART-BB medoid completion.
They must not be interpreted as an exact observed-data marginal likelihood that integrates over missing covariates.
The article discusses the inferential target, assumptions, and sensitivity analyses in detail.


## Reproducibility and limitations

This distribution contains source code, but **no NEPS individual records, prepared respondent-level data, simulation exports, execution logs, or generated results**. 
The empirical application requires authorized access to the NEPS Scientific Use Files and a valid [Data Use Agreement with LIfBi](https://www.neps-data.de/Data-Center/Data-Access/Data-Use-Agreements).
Furthermore, see [data-access instructions](docs/DATA_ACCESS.md).

The author has run the retained code in the local research environment.
Preserved results from the empirical run `sc6_20260629_103701` agree with key Bayesian findings reported in the article.
See the [execution notes](docs/EXECUTION.md).


## Quick checks (do not run a full simulation)

Run the following commands from the repository root after installing Python, R, and Julia as needed:

```powershell
# Static check (Python standard library only)
python tests/check_repository.py

# Julia dependencies and a small synthetic test
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia julia/test/conditional_medoid_smoke.jl

# R smoke tests
Rscript R/test_conditional_chib_smoke.R
Rscript R/test_ml_d2_smoke.R
Rscript R/test_ml_d3_smoke.R
```

These checks are not a substitute for the full simulation or the restricted-data analysis.
The simulation driver currently sets `RUN_MODE <- "full100"` and may launch resource-intensive computations or write to existing output directories.
Inspect its configuration and paths before running it.
For the empirical application, follow [docs/EXECUTION.md](docs/EXECUTION.md) and use only authorized NEPS data storage.


## Citation and provenance

Please cite the [published article](https://doi.org/10.1007/s10260-026-00878-5) when using the methods or code. Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

The [OSF project](https://osf.io/gjeuz/overview) originally held supplementary code from the peer-review stage.
This repository is based on the later publication working directory, not assembled from the earlier OSF files.
Historical working files remain in the author's private research archive.

*This README documents code scope and scientific provenance, not a retrospective claim that all published numbers have been independently reproduced.*
