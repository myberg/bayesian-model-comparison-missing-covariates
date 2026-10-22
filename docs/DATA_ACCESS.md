# Data access and disclosure

The empirical illustration uses **NEPS Starting Cohort 6 Scientific Use Files**.
These files are access-controlled and cannot be redistributed through this repository.
This release distributes source **code only**.
The working scripts refer to NEPS SC6 Wave-16 input version `16.0.0`.
Consult NEPS's actual data use agreement, documentation and authorized data environment before attempting a replication or publishing any aggregated output.
For details see [www.neps-data.de](https://www.neps-data.de/).
Existence of generated output in a working directory does not establish its redistribution rights.

Use environment variables (see `.env.example`):

```powershell
$env:NEPS_SC6_STATA_DIR = "<authorized-local-Stata-directory>"
$env:SC6_OUTPUT_DIR = "<private-output-directory>"
```
