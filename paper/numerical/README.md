# Numerical-study entrypoints

- `study1.R` / `study2.R`: complete fresh reproduction.
- `temp-study1-followup.R` / `temp-study2-followup.R`: the current missing-comparison runs, using existing FULL results.
- `study1.Rmd` / `study2.Rmd` and their HTML reports: notebook/report interfaces.
- `common.R`: required bootstrap; `check-*.R`: validation utilities.

Exact split-machine follow-up commands and required results are in
[publication/README.md](publication/README.md). Exploratory sampler investigations are kept outside this replication folder.

# Permanent numerical study workflows

For copy-and-run commands, result paths, completion checks and convergence interpretation, start with the [repository run guide](../../README.md#numerical-studies-mac-mini-or-linux).

`study1.R` and `study2.R` are the only experiment entrypoints. With no arguments they print plans. `--action=run` executes the complete selected study; `--action=report` rebuilds available reports without fitting. Scientific implementation lives in the installed package, not in temporary analysis scripts.

```sh
Rscript tools/setup_numerical.R
Rscript paper/numerical/study1.R --action=run --profile=paper
Rscript paper/numerical/study2.R --action=run --profile=paper
```

Study II defaults to both components. To perform only the minimal new work, use `--component=mechanism` and the selections below. The original grid retains its fixed profile; `--reps`, `--methods` and mechanism sampler controls require that explicit component. Rerunning identical commands reuses compatible outputs. The package refactor changes provenance, so pre-refactor output roots must remain separate. Pending publication displays and scientific validation remain TBD; orchestration does not manufacture missing findings.

# Numerical experiments from a fresh checkout

All numerical data, priors, seeds, fitting code and scoring code are in this repository. No Dropbox folder, saved fit, private CSV, Overleaf project, or user-specific library is required. Standard public R dependencies are restored from CRAN using `renv.lock`. This is a source-and-environment reproduction workflow; it does not claim bitwise-identical floating-point results across hardware.

## macOS or Linux setup

Install Git and R (the lock was recorded with R 4.6.0). Source-package installation may need Xcode Command Line Tools and an R-compatible Fortran compiler on macOS; Linux needs the R development toolchain, C/C++/Fortran compilers and CMake (for example `r-base-dev build-essential gfortran cmake` on Debian/Ubuntu). The dependency installer will report missing system requirements rather than substituting an implementation.

```sh
git clone https://github.com/sj156/eivGP.git
cd eivGP
Rscript tools/setup_numerical.R
Rscript tools/setup_numerical.R --check
Rscript paper/numerical/check-portability.R
Rscript paper/numerical/check-numerical.R
```

The setup restores the versioned dependency lock into `.r-library/`, then installs the checkout's `eivGP` source there. Runners locate that library and the repository from their own paths, including when invoked from another working directory. `EIVGP_R_LIBRARY` overrides the library location. `--install-local` reinstalls source using already verified dependencies; it does not download dependencies. `--check` checks dependency versions and source identity without running simulations. Keep the Git commit and any local patch with a publication run. The macOS/Linux GitHub Actions workflow runs these small validation checks; publication experiments are never launched by CI.

## Minimal design; complete reproduction

The additional design has only two settings: Study II, Gaussian four-proxy measurement, calibration 20 and 50, on 50 shared replications. Default method selection is FULL, CUT, UC-GP, LVGP, EzGP, and the three complete-case `-CC` versions. This creates 800 score jobs: 100 FULL fits, 100 CUT jobs, 300 complete-case fits, and 150 categorical-only fits whose scores are reused at both calibration levels. Each CUT job contains nested conditional GP fits, so job count is not a runtime estimate.

```sh
# Safe default: print the plan. No data or models are needed beforehand.
Rscript paper/numerical/study2.R --component=mechanism --action=plan

# Generate the original fixed-seed datasets, without fitting any model.
Rscript paper/numerical/study2.R --component=mechanism --action=data

# Small paired pilot: replications 1 and 25, calibration 20 and 50.
# 16 new CUT/CC jobs plus four matched FULL references; separate output root.
Rscript paper/numerical/study2.R --component=mechanism --action=run --profile=pilot \
  --output="$HOME/eivGP-results/study2-pilot" --cores=4

# Review study2-pilot/mechanism/report.html and stage_timing.csv first.
# After pilot review, the full fresh workflow is the two scripts above.
# Do not reuse the pilot output root for paper settings.

# Regenerate summaries only.
Rscript paper/numerical/study2.R --component=mechanism --action=report --profile=pilot \
  --output="$HOME/eivGP-results/study2-pilot"
```

Data default to `~/eivGP-results/data/study2/`; results default to `~/eivGP-results/study2/mechanism/`. Missing datasets are also generated automatically by `--action=run`. Data generation always uses the original n=100, n_test=100, q=4, m=4 design and calibration grid 0/20/50/80; the mechanism component fits only 20/50. Dataset seed is `1000000 + 10000*rep + 1`, with calibration seed one higher. Sampling seeds are recorded separately. A `--profile=quick` run keeps this DGP but uses deliberately inadequate MCMC and optimizer budgets, in a separate `--output=results/quick` directory. It tests software, not the scientific claims.

To run **only the new work** when compatible old FULL/categorical results are available:

```sh
Rscript paper/numerical/study2.R --component=mechanism --action=run --reps=1:50 --cores=4 \
  --methods=CUT,UC-GP-CC,LVGP-CC,EzGP-CC --output=results/additions
```

This is the 300-complete-case-plus-100-CUT plan, to launch only after reviewing the pilot. Add `FULL` and a compatible `--full-fit-map` to rescore 100 saved FULL references under exactly the same evaluator (500 score jobs, still 400 new comparisons). A map must include each selected replication/calibration once. Do not pair these new surface errors with historical surface-grid errors.

The pilot uses the original n=100, n_test=100 datasets and original optimizer settings: four measurement chains, 6,000 iterations / 2,000 warmup; four outer draws (one interior draw per measurement chain); four conditional GP chains per outer draw, 1,500 / 500; four retained evaluation states per conditional chain; 16 prospective latent draws per point. Replications 1 and 25 are a fixed computational convenience sample, not a random performance sample. Calibration 20 and 50 contrast weak and stronger anchoring. Shared measurement draws and common test integration draws pair CUT and all three complete-case methods. All four conditional fits retain every chain draw.

Review finite failures separately from unavailable tail ESS, chain stability for response-relevant quantities, stage wall times, and failure messages. Do not use the all-parameter flag as an automatic exclusion rule. Four outer draws cannot validate a nested posterior approximation: before the 400-job run, check additional outer draws and larger test-integration budgets on these same datasets, and compare their changes with the FULL–CUT/CC contrasts. Runtime projections must scale measurement, conditional fitting and prediction separately; the generic per-job forecast does not perform that extrapolation.

 Fresh reproduction should use the full default above so it does not depend on those old results. To refit or complete only one failed job:

```sh
Rscript paper/numerical/study2.R --component=mechanism --action=run --reps=7 --calib=20 \
  --methods=LVGP-CC --cores=4
```

## Statistical protocol and pending evidence

Categorical-only competitors fit all `(Y,X,C)` rows, independently of calibration. Complete-case competitors fit `(Y,X,U,C)` on calibrated rows, retaining C and treating U as continuous. They integrate unknown test U under a shared response-free measurement posterior and common integration draws. Neither true test U nor test Y enters prospective prediction. All methods use the same full-training definition of seen/unseen category combinations; additional point flags describe complete-case support.

CUT targets `p(U_missing,psi | C,U_observed) p(theta,f | X,Y,U)` with the same measurement and response priors as FULL. For each response-free outer draw, the standard response-GP sampler runs with all training U fixed. Outer weights remain equal; responses never reweight the measurement posterior. Defaults are four measurement chains, 20,000 iterations / 5,000 warmup, 64 outer draws, and four conditional GP chains per outer draw, 4,000 iterations / 1,000 warmup. Four states per conditional chain and 64 latent draws per test point are used for evaluation. These are starting settings: **TBD** paper-profile mixing, nested Monte Carlo sensitivity, runtime and final precision.

Evaluation distinguishes:

- Observable mean m: error against refined DGP quadrature at the same test X,C.
- Future Y: RMSE, Gaussian-mixture CRPS/NLPD, 95% coverage and width, overall and seen/unseen.
- Surface f: FULL/CUT recovery at independent test X,U, preserving response interactions; this known-input target is separate from prospective prediction.
- Missing training U: FULL/CUT errors and marginal intervals, excluding calibrated rows and retaining coordinate IDs.

For m, FULL/CUT intervals integrate the full GP covariance. CC interfaces do not provide the required covariance, so their m interval fields are unavailable; predictive intervals are not relabeled. Physical scoring checks affine-rank calibration. A fitted response f(X,U,C) is not automatically a comparable f(X,U).

`status.csv` preserves failed, pending and completed jobs and diagnostic flags. `point_scores.csv`, `replicate_scores.csv`, `summary.csv` and `paired_differences.csv` retain per-point evidence and matched replication contrasts. Parameter screens flag R-hat >1.01 or bulk/tail ESS <400; all completed scores remain visible. Replication MCSE is not within-fit error. Per-fit diagnostic CSVs and `integration_stability.csv` support review; neither a passed screen nor a split-integration check proves adequate target precision. `paired_summary.csv` gives replication MCSEs for method differences; `interaction_summary.csv` gives the paired unseen-minus-seen contrast. **TBD:** target-specific checks, budget sensitivity, and final tables and figures. The manuscript marks these unfinished analyses rather than presenting anticipated gains as findings.

## Resume, transfer and optional imports

Reissue the same command after interruption. Completed score jobs and conditional-fit checkpoints are reused. Changes to methods/replication selections or cores are allowed; changes to code, dependencies, model/sampling settings or input bytes require a new output root. Record `sessionInfo.txt`, `protocol.rds`, data hashes and `git_revision.txt` with the results. Use the same checkout and locked environment on the Mac mini or Linux machine. Do not merge outputs from differing protocols by renaming their folders.

Only one process may own an output root; `--cores` parallelizes chains within it. For scheduler jobs, give each worker a distinct output root and replication selection. An abrupt termination can leave `RUN.lock` and a job `.lock`: inspect their PID files, confirm the processes have stopped, and remove those exact locks before resuming. Caught failures release locks automatically.

Saved FULL fits are an optional optimization, not a requirement. `--full-fit-map=/absolute/path/map.csv` accepts columns `rep,calib,path` (one existing RDS per selected FULL job); data, calibration, kernel, priors, sampler and physical scaling are checked. Use absolute map paths. Supply the same map from the first launch if imported FULL and CUT/CC are to share one protocol root. Otherwise omit it and refit FULL. Rescore FULL with the new evaluator before comparing f: historical surface-grid scores use a different evaluation measure.

## Original univariate and multivariate grids

The default Study I workflow runs its original univariate grid. The default Study II workflow runs both its original multivariate grid and the targeted mechanism comparisons. Select `--component=core` or `--component=mechanism` to run or report one part independently; Study I supports core only.

```sh
Rscript paper/numerical/study1.R --action=plan --profile=paper --cores=4
Rscript paper/numerical/study2.R --action=plan --profile=paper --cores=4
# Change plan to data to generate, or to run to generate, fit and aggregate.
```

Defaults are `~/eivGP-results/data/study1`, `~/eivGP-results/data/study2`, and separate `~/eivGP-results/study1`, `~/eivGP-results/study2` run roots, with `core/` and (Study II) `mechanism/` subdirectories. An explicit `--output` sets the parent root for both components. Explicit `--data-root` and `--output` remain available. The old positional report APIs and notebooks are optional readers for historical exports; they are not prerequisites for these fresh-run commands. Local data, fit and build directories are ignored by Git. No experimental output is written to Overleaf.

## Output retention

Study fit checkpoints use lossless compact storage: pooled copies of chain
histories are reconstructed on read. The mechanism component uses the same storage.
CUT retains every conditional fit, including passing fits, together with selected prediction states and diagnostic CSVs. The selected states accelerate prediction but do not replace chain histories. Full EIV fits and shared measurement traces are retained. No existing archives are cleaned automatically.
Imputation scores preserve errors, coverage and interval widths; prediction
files omit duplicate training-U draws. Use eivGP::read_fit() for compact fits.

These changes reduce duplication, not the scientific iteration budgets. A
publication-scale disk-size bound has not been measured; full chain histories
can require substantial space. Measure a representative replication before
allocating the full run.


## Checkpoints and runtime forecasts

With no path overrides, outputs go to `~/eivGP-results/study1/` and
`~/eivGP-results/study2/`; generated data go to `~/eivGP-results/data/<study>/`.
`--output` and `--data-root` take precedence over `EIVGP_RESULTS_ROOT` and
`EIVGP_DATA_ROOT`, which take precedence over these defaults. Existing results
are not moved. To resume an older run, explicitly select its original output
and data directories and use the compatible code/package version.

Both study entrypoints save completed fits immediately and reuse them on rerun.
The default partial-checkpoint schedule is **every three hours, enabled after six
hours of actual script runtime**. The first eligible save is at the first iteration
boundary after hour six, followed by script-clock boundaries at hours nine,
twelve, and so on. Completed chains also save once this six-hour gate is open.
The gate uses elapsed runtime, not a prediction of total duration, and resets on
each command invocation. Before it opens, interrupted fits restart unless a
compatible checkpoint from an earlier invocation exists. Existing compatible
checkpoints can still be loaded. Each completed CUT conditional fit is cached
separately. Timing logs and forecasts still update about every minute.

Defaults are `--checkpoint-seconds=10800 --checkpoint-after-seconds=21600
--checkpoint-iterations=Inf`. The iteration-count trigger is disabled. An explicit
`--checkpoint-seconds` or `--checkpoint-iterations` override defaults the delay to
zero unless `--checkpoint-after-seconds` is also supplied. Set both intervals to
`Inf` for completed-fit-only operation: no partial checkpoints are written or
loaded, including at individual chain completion. Direct package sampler calls
also use completed-fit-only behavior unless checkpoint interval options are set.
Checkpoint saves occur at iteration boundaries, so a slow iteration can delay a
save. Third-party optimizers and prediction calculations retain their existing
completed-stage caching. Existing checkpoint files are never deleted merely
because the schedule changes.

The following checkpoint details apply only when partial-fit recovery is enabled.
`<output>/checkpoints/<fit-identity>/chain-XX.rds` stores all generated draws, current
state, RNG state, and (for GP chains) computational cache and counters. Inputs,
sampling settings, seeds, package module hashes, R version/platform and RNG kind
identify a fit. Reissuing the same command restores saved iterations and reuses
completed chains. Writes use a temporary file and same-directory rename, leaving
the previous checkpoint intact until replacement succeeds. This protects against
partial writes; it is not a guarantee against filesystem or disk failure.
There is one rolling checkpoint per chain, not one file per save. Checkpoints
use lossless gzip compression and contain only generated draw/trace rows, omitting
unused allocation for future iterations. Restart restores those rows into the
full in-memory arrays. Every generated draw is retained; no thinning is applied.
A replacement temporarily requires space for both the old and new compressed file.

When a full GP fit is published through the numerical artifact writers, the
writer reads it back and checks exact equality, all chain draws, model inputs and
settings, initial/terminal states, RNG, computational cache and sampler counters. Only then does it
remove that fit's redundant `chain-XX.rds` files under the checkpoint lock.
The identity file is retained for provenance. The full fit (including all conditional CUT chains, diagnostics
and continuation state), timing histories, and a small `retired.rds` receipt
remain. Partial fits, continuation-segment checkpoints, selected prediction states, scores, failed read-back checks,
locked directories and checkpoints outside the configured root do not trigger
retirement. If no full fit has been published, restart files remain necessary.
Study II measurement checkpoints are compressed but retained: its analysis
artifact does not preserve the full terminal sampler state. Cleanup adds a
read-back verification cost at fit publication; checkpoint compression also costs
CPU time, in exchange for reduced storage and I/O.

This policy applies to newly generated, linked artifacts. It does not scan or
delete older output trees, infer completion from a filename, or migrate historical
runs. Existing portable-drive results need an inventory and verification against
their final fits before cleanup. The reader supports older full-allocation
checkpoint layouts, but the existing code/protocol compatibility checks still
apply: do not switch a running experiment to new code or bypass its identity
checks merely to shrink its checkpoints.

After a crash, a `<fit-identity>/LOCK/owner.csv` may remain. Verify its recorded
host and PID are no longer running before removing that exact `LOCK` directory.
Also check mechanism `RUN.lock` and job `.lock` directories as described above.
Never run two processes against the same output root. Do not change code,
dependencies, inputs or scientific settings to resume a run. Checkpoint cadence
and core count may change. The new implementation does not reconstruct partial
chains from older runs that never saved them.

Each chain owns `chain-XX-timing.csv` (latest) and
`chain-XX-timing-history.csv` (history). They record UTC update time, current and
saved iterations, warmup/sampling/completion phase, elapsed active wall seconds,
remaining seconds and estimated finish time. Updates occur after 60 seconds or
a checkpoint, at iteration boundaries. Elapsed time carries forward from the
last saved checkpoint on restart; downtime and lost unsaved work are excluded.
The chain forecast extrapolates its observed average iteration cost, so early
warmup estimates are provisional and do not include queued chains, diagnostics,
prediction or other fits. These are runtime estimates, not convergence forecasts.

Core replication timing lives beside each replication status CSV in its
`.timing/` directory; mechanism job timing lives in `mechanism/timing/`.
Per-task RDS records preserve completed-attempt durations and accumulated
recorded task seconds, including caught failures. `forecast-<PID>.csv` and its
history update on task starts/completions/failures/cache hits. Use the latest UTC
update across these worker-owned files. Core forecasts cover that replication
batch only, using mean completed-replication wall time divided by configured
replication workers; mechanism forecasts use separate mean times by method and
serial job scheduling. Running tasks are conservatively counted as full remaining
tasks. A forecast remains NA until every unfinished group has an observed completed
duration. Cache-load time never substitutes for model runtime. Predictions are
rough planning estimates: heterogeneous tasks, restart work, scheduling and shared
measurement costs can change them. No single batch forecast claims to cover all
later study cells or competitors. Completed final batch forecasts report zero
remaining work; task seconds are summed work, not elapsed end-to-end wall time
when workers run in parallel.

```sh
# Discover current chain and batch forecasts under the default roots.
find "$HOME/eivGP-results" -type f \
  \( -name 'chain-*-timing.csv' -o -name 'forecast-[0-9]*.csv' \)
```

## Notebook parity

The `.Rmd` files can read historical merged exports when `component` is empty, or call the same package workflow as the two batch scripts when it is set. For a saved Study II pilot, render `study2.Rmd` with `component="mechanism"`, `profile="pilot"`, `action="report"`, `results_root` equal to the batch `--output`, and `options=list(`data-root`=..., `full-fit-map`=...)` matching the original launch. Omit the FULL map if none was used. Rendering reports never silently changes action to run. `study1.R` includes its original grid and complete-case published predictors; `study2.R` includes its core grid and the CUT/complete-case component. Both are included automatically with `--component=all` (the default). No third simulation script is required.


## Two versions of the published competitors in both studies

The original UC-GP/LVGP/EzGP fits use every `(Y,X,C)` record, ignore calibration entirely, and are reused across calibration sizes. Their `-CC` versions use only calibrated `(Y,X,U,C)` records, with `U` appended to the continuous inputs and `C` retained as categorical. Unknown test inputs are integrated over common response-free measurement draws; no test `Y` or true test `U` enters prediction. The latter integration is an adaptation, not a claim about the original published algorithms. All three retain their author-package fits and optimizer status. Missing category support and numerical failures remain explicit; no failed fit is replaced by an original-method result.

Study I now has `--component=complete-case`. Its paper grid adds **600 CC fits** (two regimes, 50 replications, calibration 20/50, three methods), separate from the **400 new Study II jobs** (100 CUT and 300 CC). `--component=all` runs the core original-method grid plus these additions through each study's single script. Use `--component=core` to reproduce only the prior core grid. No paper batch is launched automatically after a pilot.

```sh
# First inspect and then run 12 paired pilot jobs: both regimes, replication 1,
# calibration 50, three original and three complete-case predictors.
Rscript paper/numerical/study1.R --component=complete-case --profile=pilot --action=plan
Rscript paper/numerical/study1.R --component=complete-case --profile=pilot --action=run --cores=4 --output="$HOME/eivGP-results/study1-cc-pilot"
# Saved results only; no refitting. Use the same profile and scientific settings.
Rscript paper/numerical/study1.R --component=complete-case --profile=pilot --action=report --output="$HOME/eivGP-results/study1-cc-pilot"
# After review, run only the new Study I fits (originals remain in core).
Rscript paper/numerical/study1.R --component=complete-case --profile=paper --action=run --output="$HOME/eivGP-results/study1-cc-paper"
```

Study I accepts `--eta=0,1`, `--reps`, `--calib=20,50`, and `--methods` within this component. The pilot uses four response-free chains with 6000 iterations/2000 warmup and four balanced outer states with 16 prospective latent draws each; paper starts with 20000/5000, 64 states and 64 latent draws. Optimizer controls remain the paper controls in the pilot; `quick` is only a software smoke check. The response-free deterministic threshold sampler uses `U~N(0,1)` and ordered uniform thresholds on `[-8,8]`, as in the Study I diagnostic baselines; it is **not** the Study II probit model or a matched Study I FULL-CUT experiment. It saves complete measurement chains at the end of the stage; unlike the Study II CUT sampler it has no within-chain restart checkpoint. Completed measurement and native GP stages are cached for retry.

Per-regime `status.csv`, `summary.csv`, `paired_summary.csv`, `report.html`, native fits, measurement chains, diagnostic CSVs and stage times are saved under `complete-case/eta0_balanced/` and `eta1_balanced/`. Both regimes are reported separately. A one-replication pilot has no replication MCSE and cannot establish comparative performance. Review mixing and latent-integration sensitivity before increasing the grid. The full `m` point estimate uses all latent components; response scores by default use four evenly spaced latent components per outer state (configurable with `--cc-response-latent`), a budget that must also be checked. Complete-case mean credible bands are unavailable because native adapters do not provide the integrated joint covariance.

For the Study I notebook use `component="complete-case"`, `profile="pilot"`, `action="report"`, the batch output root, and any matching selections in `options`. PI-GP/CC-GP remain supplementary diagnostics; neither is required as a primary benchmark once the original and complete-case versions are reported, and neither replaces the matched FULL-CUT feedback comparison.


## Stronger paired pilot and external output storage

Keep `--profile=pilot` and select only the paired pilot replications. No pilot command switches to a full grid automatically. Large outputs and `checkpoints/` follow `--output`; frozen data follow `--data-root`. On an external drive, verify that it is mounted before launching and use a new output root for changed scientific settings. Preserve the initial pilot outputs for comparison.

`--latent` controls integration for the observable mean. `--gp-response-latent` controls prospective U draws per retained GP state for FULL/CUT response scores (default 1); `--cc-response-latent` controls draws per outer measurement state for CC response scores (default 4). Both must be positive and no larger than `--latent`. All components retain equal weights. Increasing either response budget changes no fitted posterior and does not create mean credible intervals for the native competitors. The settings and package identity are recorded and must match on resume/report.

```sh
# Set PILOT_STORE to a mounted local/external results directory outside Overleaf.
# Study I: both regimes, matched replications and calibration levels; 48 score
# jobs, with original-method fits cached across the two calibration levels.
Rscript paper/numerical/study1.R --component=complete-case --profile=pilot --action=run \
  --eta=0,1 --reps=1,25 --calib=20,50 --cores=2 \
  --iter=12000 --warmup=4000 --outer=16 --latent=128 --cc-response-latent=16 \
  --output="$PILOT_STORE/study1" --data-root="$PILOT_STORE/data/study1"
# Study II: only 20 score jobs, including four FULL references; no core grid.
# Supply --full-fit-map=/absolute/map.csv to rescore compatible saved FULL fits.
Rscript paper/numerical/study2.R --component=mechanism --profile=pilot --action=run \
  --reps=1,25 --calib=20,50 --cores=4 --iter=100000 --warmup=20000 \
  --outer=16 --inner-iter=4000 --inner-warmup=1000 --inner-draws=4 \
  --latent=128 --gp-response-latent=4 --cc-response-latent=16 \
  --output="$PILOT_STORE/study2" --data-root="$PILOT_STORE/data/study2"
```

These are numerical-accuracy pilots, not performance estimates or validated final budgets. The longer Study II measurement budget targets the low ESS found in the initial pilot. Review target/measurement diagnostics, compare nested outer-state subsets, and repeat prediction with independent integration seeds before setting production budgets. Retain finite failures and unavailable diagnostics separately. The Study I deterministic measurement traces are shared in law across its two response regimes; identical diagnostics are not independent evidence.

### Default mixing figures

Both `study1.R` and `study2.R` (and their notebooks) save three small ggplot
figures after `--action=report` or `--action=run`: chain agreement versus
bulk ESS, MCSE/posterior SD, and unavailable diagnostic counts. Each selected
component writes a separate `mixing-plots/` folder containing PDF and PNG
versions plus a diagnostic-source manifest. Set `--mixing-plots=false` to skip
them. Planning and data preparation never fit or plot chains. Plotting reads
saved parameter-level CSVs; historical exports containing only gate summaries
are explicitly reported as insufficient rather than inventing diagnostics.

Reusable package APIs:

```r
eivGP::plot_mixing(diagnostics, type = "overview") # also precision / coverage
eivGP::save_mixing_plots(results_directory)
eivGP::plot_mcmc(fit, parameters = c("V", "r"), type = "trace")
# also rank / acf / running_mean; supply named chain-vector lists for targets
```

These functions return ordinary ggplots (the save helper returns file paths).
Trace display subsampling never changes diagnostic calculations or stored draws.
Scientific-target integrations are not recomputed for plots. Keep raw-coordinate,
scientific-target, and pilot results distinct; passing reference lines establishes
neither convergence nor integration accuracy. Reinstall the package after updating
this checkout; a running pilot's frozen library should remain unchanged.

## Publication completion (2026-10-19)

Both study entrypoints support `--component=publication --stage=pilot|integration|full`.
See [portable publication workflow](publication/README.md). Default action is plan;
no full simulation has been launched by this update.

### Shared CPU budget

`--cores` supplies one budget for the entire numerical workflow, including
replication workers, the four chains per fit, posterior prediction, diagnostic
summaries, scoring and mixing figures. Waiting coordinators release their slots.
Independent replication groups keep ownership of their shared measurement and
prospective-input caches; no method races another method to write those caches.
Set BLAS/OpenMP thread counts to one before R starts, as described in the root
README. Windows retains serial execution.

Prediction stream schema 1 uses fixed blocks of eight selected posterior states
and seeds independent of core count. Serial and parallel results agree within
this schema; older serial prediction realizations may differ. Existing result
and checkpoint identities remain enforced. Use a new output root when the
recorded identity is incompatible; do not delete provenance files to bypass it.


## Verified cleanup of existing Study II CUT checkpoints

Use the package function `prune_study_checkpoints()` on the directory containing
both `checkpoints/` and `mechanism/` (for follow-up runs, the `reference/` folder).
It audits old checkpoint trees that predate automatic retirement. The default is
a dry run: it writes an audit only, without changing the experiment output tree.
Choose a fresh local audit directory outside the experiment tree on every call.

```r
library(eivGP)
root <- "/path/to/study2-followup/reference"
prune_study_checkpoints(root, audit_dir = "~/checkpoint-audits/dry-run-01")
prune_study_checkpoints(root, dry_run = FALSE,
                       audit_dir = "~/checkpoint-audits/cleanup-01")
```

The destructive call repeats verification; it does not blindly trust a previous
manifest. A seed or matching filename only identifies a candidate. The verifier
compares the data, normalized kernel, priors, sampler controls, budgets, every
saved draw in all chains, initial states, terminal states, RNG, cache and counters.
Squared-exponential kernels ignore their unused Matérn smoothness argument. A
unique exact match is required. File hashes are rechecked immediately before
removal, while holding the mechanism run lock and the individual sampler lock.
Existing locks are respected, never broken automatically.

Only verified `chain-XX.rds` files are removed. All full conditional fits,
measurement checkpoints, incomplete/mismatched/ambiguous work, continuation
segments, identities, timing logs, predictions and scores remain. Current automatic
retirement uses the same verifier. This legacy cleanup targets the Study II CUT
layout; other checkpoint types are retained. Code/protocol compatibility rules for
resuming experiments still apply; running this offline utility does not upgrade,
rewrite or relabel their saved scientific results.

`manifest.csv` records each directory's status, preservation reason, matched fit,
verified bytes and removed bytes. `receipts/*.rds` record exact filenames, hashes,
original sizes and elapsed chain times before deletion and removal status afterwards.
If interrupted, inspect those receipts and retry with a fresh audit directory.
Never manually remove a run/sampler lock until all owning processes have stopped.
