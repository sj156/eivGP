# Numerical simulations

Study-specific scripts retain their Study I / Study II names. `setup_study*`
configures the low-level drivers; `00_synthetic_data.R` handles frozen inputs;
`02_*` fits replicated datasets; `04_*` implements ablations. The `01_*` scripts
produce representative figures. Publication recovery and competitor caching
are repository services, not installed package functions.

The shared simulation implementation is split under `helpers/`:

- `design.R`: configuration validation, resource allocation, study cells and constructors.
- `data.R`: task plans, preflight, manifests, source hashes, frozen-data checks and generation.
- `execution.R`: per-cell controls and replicated study execution.
- `gates.R`: diagnostic and output completeness checks.
- `summaries.R`: aggregation, paired comparisons and publication summaries.
- `workflow.R`: top-level simulation orchestration and dataset identity.

Load all helpers through `source("codes/simulation_helpers.R")` from the
repository root. CLI launchers in `codes/cli/` are the normal entry points.
