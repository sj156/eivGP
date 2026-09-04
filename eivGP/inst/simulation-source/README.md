# Mixed-input EIV-GP computational workflow

For the publication simulations, start with
[`PUBLICATION_SIMULATIONS.md`](PUBLICATION_SIMULATIONS.md). The current
publication interface is one helper bundle (`reproduction_workflows.R`) and one
master script per study (`run_study1_simulation.R` and
`run_study2_simulation.R`). Both masters default to a read-only dry run.

The posterior engines are shared with the planned R package, while publication
orchestration is intentionally narrower. Run the following from this directory:

```sh
Rscript run_study1_simulation.R
Rscript run_study2_simulation.R
```

Those commands are dry runs. Generate frozen data, smoke-test the complete
workflow, and launch the archival runs by setting `MIXEDGP_RUN_MODE` to
`data`, `smoke`, and `publication`, respectively. For example:

```sh
MIXEDGP_RUN_MODE=publication MIXEDGP_WORKERS=8 \
  Rscript run_study1_simulation.R
```

Set `MIXEDGP_R_LIBRARY=/absolute/path/to/R-library` when the audited package
library is not stored beside the project root. Optional `MIXEDGP_DATA_ROOT`
and `MIXEDGP_OUTPUT_ROOT` variables place large workstation artifacts outside
the synchronized Overleaf tree without changing the masters.

`reproduction_workflows.R` is the single source of truth for design cells, seeds,
task eligibility, frozen-data checks, dispatch, gates, aggregation, and paired
comparisons. The masters contain only the user-editable run configuration and
one function call. The scientific posterior engines remain separate because
they are also the package source of truth; they are not duplicated in the
simulation layer.

## Canonical modules

- `model_api.R`: package-facing EIV-GP and competitor calls.
- `reproduction_workflows.R`: publication design and orchestration bundle.
- `run_study1_simulation.R` and `run_study2_simulation.R`: publication masters.
- `reproduction_compat.R`: package-facing and historical unified workflow.
- `core_parallel.R`: deterministic `mclapply`/`mcmapply` wrappers.
- `core_numerics.R`: shared linear algebra, diagnostics, and scoring utilities.
- `reproduction_data.R`: versioned frozen datasets and checksum manifests.
- `model_univariate.R` and `model_multivariate.R`: reusable posterior engines.
- `competitors.R`: audited public-package adapters shared by both studies.
- `08_published_competitor_validation.R`: exact small-data fit and prediction
  checks for every published competitor adapter.
- `09_experiment_design_validation.R`: exact paired-design and nested-proxy
  invariance checks shared by both studies.

The `run_study1_all.R` and `run_study2_all.R` files are compatibility aliases
for the audited masters. `run_mixedgp_experiment()` remains available for
package examples and historical workflows, but it is not the publication
simulation entry point.

## Reproducibility gates

The study runners now enforce the formerly manual validation sequence. Before
every smoke or publication fit they run the exact checks in
`08_published_competitor_validation.R` and
`09_experiment_design_validation.R`, archive the results under the
fingerprinted run's `config/` directory, and stop on a failed invariant.
Publication mode additionally requires all prespecified competitors to pass
their small real-fit validation.

The end-to-end workflow is therefore:

1. run each master without environment overrides and inspect the dry-run
   specification and package preflight;
2. run both masters with `MIXEDGP_RUN_MODE=smoke`; the validators execute
   automatically before the non-reportable fits;
3. run both masters with `MIXEDGP_RUN_MODE=publication`; the validators execute
   automatically and publication mode activates strict competitor and MCMC
   gates;
4. archive manifests, validation tables, raw RDS results, diagnostics, package
   versions, and session information together.

Generated datasets live under `../data-synthetic/publication-v2/`. Their
CSV/RDS manifests record schema versions, seeds, sizes, and MD5 checksums.
Before fitting, the masters also verify the common-random-number and nested-
proxy contracts. Fitting stops when a required artifact is absent or altered;
it never silently regenerates data in a fit-only run.
