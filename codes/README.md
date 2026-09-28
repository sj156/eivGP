# Supporting R code

`core/` is the canonical source used to generate the installable `eivGP/`
package through `litr/`. Other collections implement repository workflows.
The model posterior, sampler, study settings, and fitting budgets are not
changed by this directory reorganization.

| Collection | Responsibility |
| --- | --- |
| `core/` | Model utilities, current sampler, public API, MCMC diagnostics and published competitors |
| `simulations/` | Synthetic data, study drivers, setup, ablations, caching and publication recovery |
| `reporting/` | Reporting saved fits and combining experiment outputs |
| `applications/` | ADNI and ocean analysis support, shared paths, dependency setup, and runners |
| `cli/` | Repository-root command-line entry points for numerical workflows |
| `validation/` | Computational checks and design pilots |
| `tests/` | Automated regression tests |
| `setup/` | Project library and thread configuration |
| `legacy/` | Historical workflows and case-study artifacts retained for provenance |

Two stable source entry points remain here:

```r
source("codes/load_mixedgp.R")       # Canonical model modules in dependency order
source("codes/simulation_helpers.R") # Study design, data, execution and summaries
```

`simulation_helpers.R` loads `simulations/helpers/` into the caller's environment.
The helper files separate design, frozen-data management, execution, diagnostic
gates, summaries, and workflow orchestration. Edit the appropriate helper rather
than adding another copy. Source the loader, not an individual helper file.

The historical `00_study1_functions.R` and `00_study2_functions.R` names in
`core/` describe their origins; both are shared model implementations. Load the
current sampler and API after these modules using `load_mixedgp.R`. Application
support uses these same canonical files; duplicate model copies were removed.

Use launchers from the repository root, for example:

```sh
EIVGP_RUN_MODE=dry_run Rscript codes/cli/run_publication_study.R study1
MIXEDGP_CORE_BUDGET=4 Rscript codes/cli/run_development_study.R study2 plan
Rscript codes/applications/run_real_data.R adni --cores 12 --plan
Rscript codes/tests/test_v031_updates.R
```

Low-level study and pilot scripts that use relative `source()` calls retain their
`codes/` working-directory convention: from `codes/`, source
`simulations/...` or `validation/...`. Prefer the documented launchers for full
workflows. `run_study1_all.R` and `run_study2_all.R` remain compatibility aliases
inside `simulations/`.

All active analysis notebooks are in `../replication/`; package-building
notebooks remain in `../litr/`. Generated results stay in the configured output
directories, usually `../reproduction/` or `../real_data_application_outputs/`.
Historical material in `legacy/` is not a supported current analysis entry point.

For details, see [publication simulations](PUBLICATION_SIMULATIONS.md),
[Study II](README_study2.md), and
[applications](applications/REAL_DATA_APPLICATIONS.md).
