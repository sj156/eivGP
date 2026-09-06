# Publication simulation workflow

September 5 computational status: publication release remains on hold pending
adequately mixed posterior runs. Both drivers now require target diagnostics
as well as the applicable raw/invariant checks (four chains, Rhat <= 1.01,
bulk/tail ESS >= 400). Initial chain lengths are not convergence guarantees.
Run `16_sampler_efficiency_benchmark.R` for bounded same-target computation
checks. Its short runs and saved timing tables are explicitly non-reportable.
The command examples below describe the workflow, not approval of the current
pilot results.

The publication workflow has three entry files:

- `simulation_helpers.R`: the only simulation helper/orchestration bundle;
- `run_study1_simulation.R`: the Study I master;
- `run_study2_simulation.R`: the Study II master.

The posterior engines remain in `00_study1_functions.R` and
`00_study2_functions.R` because they are also the package source of truth.
Master scripts do not source the historical representative, pilot, or unified
runner entry points.

## Safe use

Run these commands from `revision/codes`. Both masters default to a read-only
dry run:

```sh
Rscript run_study1_simulation.R
Rscript run_study2_simulation.R
```

Generate and checksum the frozen synthetic datasets without fitting:

```sh
MIXEDGP_RUN_MODE=data MIXEDGP_WORKERS=8 Rscript run_study1_simulation.R
MIXEDGP_RUN_MODE=data MIXEDGP_WORKERS=8 Rscript run_study2_simulation.R
```

Run a non-reportable workflow smoke test:

```sh
MIXEDGP_RUN_MODE=smoke MIXEDGP_WORKERS=4 Rscript run_study1_simulation.R
MIXEDGP_RUN_MODE=smoke MIXEDGP_WORKERS=4 Rscript run_study2_simulation.R
```

Run the frozen publication design on the workstation:

```sh
MIXEDGP_RUN_MODE=publication MIXEDGP_WORKERS=8 Rscript run_study1_simulation.R
MIXEDGP_RUN_MODE=publication MIXEDGP_WORKERS=8 Rscript run_study2_simulation.R
```

`MIXEDGP_STAGES=data`, `fit`, or `fit,aggregate` may be used to resume a
specific stage, but only together with `MIXEDGP_RUN_MODE=smoke` or
`MIXEDGP_RUN_MODE=publication`. For example, resume a reportable fit and
aggregation with:

```sh
MIXEDGP_RUN_MODE=publication MIXEDGP_STAGES=fit,aggregate \
  MIXEDGP_WORKERS=8 Rscript run_study2_simulation.R
```

A fit-only run verifies every selected manifest entry, generator version,
design field, and MD5 checksum before dispatch. Replications are parallel;
posterior chains are kept serial within each worker to avoid nested
parallelism. BLAS/OpenMP threads are set to one during a run.

## Frozen designs

Study I is a paired mechanism study. The balanced primary cells use
heterogeneity strengths `eta = 0, 0.5, 1` with common random numbers and the
same observable mean `m(x,c)`. The `eta = 0` cell is a category-sufficient
negative control, not a mathematical special case of LVGP or EzGP: its latent
surface is piecewise constant in `u`. An imbalanced-threshold `eta = 1` cell is
an appendix sensitivity.

Study II tests the shared-state advantage. The primary proxy-dimension cells
use nested `q = 2, 3, 4` measurements; `q = 4` also carries the calibration
curve. Additive-response and high-uncertainty controls are focused at 10%
calibration. Logistic score errors are a misspecification sensitivity, and
`q = 6` is a sparse-pattern stress test. PI-GP and CC-GP are run at the focused
calibration size in both the additive and high-uncertainty controls, making the
integration-versus-plug-in claim directly testable rather than relying only on
categorical competitors.

## Fair method-by-target comparisons

- `m(x,c)` and prediction of `Y* | x*, c*`: EIV-GP, UC-GP, LVGP, and EzGP.
  Published competitors' direct package predictive means are used for
  `m(x,c)`; noisy predictive draws are used only for response prediction.
- `f(x,u)`: EIV-GP, PI-GP, CC-GP, and the infeasible Full-U GP benchmark.
  Categorical mixed-input competitors are not assigned this target.
- latent `U`: EIV-GP and the response-free ordinal measurement model, for both
  missing training values and prospective subjects where implemented.
  Physical-coordinate scores are deliberately not reported without anchoring
  calibration.

Every run writes its resolved configuration, task plan, estimand-method
matrix, package preflight, source checksums, frozen-data manifest, cell status,
diagnostic gates, raw outputs, result checksums, and session information under
a fingerprinted, resumable run directory. It also writes paired EIV-versus-
competitor differences for every eligible task, Monte Carlo standard errors,
and the prespecified mechanism/design contrasts. Positive reported advantages
always favor EIV-GP. Publication mode stops when a
required package, fit, convergence gate, checksum, or failure-rate gate fails.

The key combined files are:

- `method_comparisons_paired_raw.csv` and
  `method_comparisons_paired_summary.csv` for prediction, `m(x,c)`, `f(x,u)`,
  and latent-`U` targets;
- `design_contrasts_paired_raw.csv` and
  `design_contrasts_paired_summary.csv` for the Study I heterogeneity contrast
  and the Study II interaction, uncertainty, and proxy-dimension contrasts;
- `all_raw_outputs.rds` and `result_manifest.csv` for lossless downstream
  reporting and checksum verification; non-tabular per-cell metadata are
  preserved separately in `cell_non_tabular_outputs.rds`.

The dry run checks all runtime packages before computation: `ggplot2`,
`dplyr`, `tidyr`, `knitr`, and `posterior` for both studies;
`TruncatedNormal` for Study II; and the published-competitor packages `kergp`,
`LVGP`, and `EzGP`. Source-file hashes and all of these package versions enter
the run fingerprint, so changed code or software cannot silently reuse an old
publication cache.

If the dependencies are kept in a project library outside the Overleaf tree,
point the masters to it without editing either script:

```sh
MIXEDGP_R_LIBRARY=/absolute/path/to/R-library \
  Rscript run_study1_simulation.R
```

For long runs, `MIXEDGP_DATA_ROOT` and `MIXEDGP_OUTPUT_ROOT` may point to a
local high-throughput disk; both resolved absolute locations are recorded in
the run metadata. Leaving them unset keeps data and results under the project.
