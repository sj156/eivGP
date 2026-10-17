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

# First paper-profile replication: assess cost, convergence and integration error.
Rscript paper/numerical/study2.R --component=mechanism --action=run --reps=1 --cores=4

# Then complete the same protocol. Finished jobs are skipped automatically.
Rscript paper/numerical/study2.R --component=mechanism --action=run --reps=1:50 --cores=4

# Regenerate summaries only.
Rscript paper/numerical/study2.R --component=mechanism --action=report
```

Data default to `~/eivGP-results/data/study2/`; results default to `~/eivGP-results/study2/mechanism/`. Missing datasets are also generated automatically by `--action=run`. Data generation always uses the original n=100, n_test=100, q=4, m=4 design and calibration grid 0/20/50/80; the mechanism component fits only 20/50. Dataset seed is `1000000 + 10000*rep + 1`, with calibration seed one higher. Sampling seeds are recorded separately. A `--profile=quick` run keeps this DGP but uses deliberately inadequate MCMC and optimizer budgets, in a separate `--output=results/quick` directory. It tests software, not the scientific claims.

To run **only the new work** when compatible old FULL/categorical results are available:

```sh
Rscript paper/numerical/study2.R --component=mechanism --action=run --reps=1:50 --cores=4 \
  --methods=CUT,UC-GP-CC,LVGP-CC,EzGP-CC --output=results/additions
```

This is the 300-complete-case-plus-100-CUT plan. Fresh reproduction should use the full default above so it does not depend on those old results. To refit or complete only one failed job:

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

`status.csv` preserves failed, pending and completed jobs and diagnostic flags. `point_scores.csv`, `replicate_scores.csv`, `summary.csv` and `paired_differences.csv` retain per-point evidence and matched replication contrasts. Parameter screens flag R-hat >1.01 or bulk/tail ESS <400; all completed scores remain visible. Replication MCSE is not within-fit error. Per-fit diagnostic CSVs and `integration_stability.csv` support review; neither a passed screen nor a split-integration check proves adequate target precision. **TBD:** target-specific checks, budget sensitivity, final tables and figures, and the unseen-minus-seen contrast with paired uncertainty. The manuscript marks these unfinished analyses rather than presenting anticipated gains as findings.

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
CUT stores selected prediction states and diagnostic CSVs for each conditional
fit, then removes that fit only when its MCMC screening passes. Failed fits
remain available for investigation. Selected states support deterministic
prediction reconstruction, but cannot resume the original conditional chains;
rerunning those chains requires a fresh output root. No existing archives are
cleaned automatically. Full EIV fits and shared measurement traces are retained.
Imputation scores preserve errors, coverage and interval widths; prediction
files omit duplicate training-U draws. Use eivGP::read_fit() for compact fits.

These changes reduce duplication, not the scientific iteration budgets. A
publication-scale disk-size bound has not been measured; failed diagnostics can
still retain substantial histories. Measure a representative replication before
allocating the full run.


## Checkpoints and runtime forecasts

With no path overrides, outputs go to `~/eivGP-results/study1/` and
`~/eivGP-results/study2/`; generated data go to `~/eivGP-results/data/<study>/`.
`--output` and `--data-root` take precedence over `EIVGP_RESULTS_ROOT` and
`EIVGP_DATA_ROOT`, which take precedence over these defaults. Existing results
are not moved. To resume an older run, explicitly select its original output
and data directories and use the compatible code/package version.

Both study entrypoints enable disk checkpoints for the shared GP sampler and
Study II CUT measurement sampler, including warmup. Every chain saves at the
first iteration boundary after 300 seconds or 1,000 iterations since its last
save, and at completion. Override with `--checkpoint-seconds=120` and/or
`--checkpoint-iterations=500`. A single slow iteration can exceed the requested
time interval. Third-party optimizers and prediction calculations retain their
existing completed-stage caching; they do not gain within-iteration restart.

`<output>/checkpoints/<fit-identity>/chain-XX.rds` stores all draws, current
state, RNG state, and (for GP chains) computational cache and counters. Inputs,
sampling settings, seeds, package module hashes, R version/platform and RNG kind
identify a fit. Reissuing the same command restores saved iterations and reuses
completed chains. Writes use a temporary file and same-directory rename, leaving
the previous checkpoint intact until replacement succeeds. This protects against
partial writes; it is not a guarantee against filesystem or disk failure.
Keep sufficient free disk space: raw restart files retain full chain arrays,
in addition to final analysis artifacts, and a write temporarily needs room for
both versions. No draws are thinned for checkpointing.

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
