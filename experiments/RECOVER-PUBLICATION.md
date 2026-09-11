# Recover the original 50-dataset publication comparisons

Run this from the repository containing the revised `codes/`, `eivGP/`, and `experiments/` directories. Run **Study II on Linux and Study I on the Mac mini**, each with its own original frozen `.rds` inputs and publication archive. Do not regenerate datasets to replace failures. Recovery does not modify the archives, original reports, or the Overleaf paper.

## Commands

From the updated repository root on **Linux (Study II)**, with the usual `reproduction/` layout:

```sh
Rscript --vanilla experiments/recover_publication.R plan --studies=study2
Rscript --vanilla experiments/recover_publication.R run --studies=study2 --cores=4 > recovery-study2.log 2>&1
Rscript --vanilla experiments/recover_publication.R check --studies=study2
```

On the **Mac mini (Study I)**:

```sh
Rscript --vanilla experiments/recover_publication.R plan --studies=study1
Rscript --vanilla experiments/recover_publication.R run --studies=study1 --cores=4 > recovery-study1.log 2>&1
Rscript --vanilla experiments/recover_publication.R check --studies=study1
```

The plan prints the missing IDs without fitting. Run the same `run` command again to resume interrupted work. `check` independently checks input provenance, original replication IDs, duplicates, and finite predictive scores. `run`/`check` return a nonzero exit status if any requested predictive comparison remains incomplete.

To use the SANDISK archives and an explicit frozen-input location:

```sh
Rscript --vanilla experiments/recover_publication.R plan \
  --archive-root="/Volumes/SANDISK USB" \
  --data-root="/path/to/original/reproduction/data/synthetic" \
  --output="/path/to/recovered-publication"
```

Replace `plan` with `run` or `check`, preserving the same options. The Linux mount path will differ from the macOS `/Volumes/...` path. `--data-root` is the parent of the **study1/** and **study2/** input directories. For cells needing MCMC, also retain the saved replication-1/calibration-50 fit under the source archive's `cells/<cell>/results/` tree; the supplied SANDISK archive contains these. Its priors, kernel, and sampler controls are passed explicitly to new fits, with a source checksum recorded, so changed package defaults cannot change the recovered model. The source run's checksums must match; supplying dataset identity CSVs alone is insufficient. The script neither copies nor regenerates frozen inputs.

`--studies=study2` limits recovery to Study II. Use a separate `--output` folder if changing the study selection. The default recovers both studies, including Study I EzGP failures. Study II always retains the original calibration grid, so it also repairs the calibration curves, not just calibration 50 in Table 2.

### Automatic use of the core allocation

Set `--cores=N` to the CPU allocation for this job (default 4). The runner automatically divides it across independent EIV–GP datasets and their archived MCMC chains. Competitor concurrency also defaults to the full allocation. For the 56-core Linux machine:

```sh
Rscript --vanilla experiments/recover_publication.R plan --studies=study2 --cores=56
Rscript --vanilla experiments/recover_publication.R run --studies=study2 --cores=56 > recovery-study2.log 2>&1
Rscript --vanilla experiments/recover_publication.R check --studies=study2 --cores=56
```

With four chains per dataset, this permits **14 concurrent datasets × 4 chains = 56 chain workers**, or up to 56 independent competitor fits in the separate competitor phase. The default four-core run remains one dataset × four chains. No extra chains or MCMC iterations are introduced. Dataset slots refill as soon as a task finishes; the coordinator immediately merges and saves each completed result. Each dataset retains its own fit and predictive checkpoints.

Optional limits:

- `--chain-workers=4` limits concurrent chains within a dataset; the original four chains still run when fewer workers are requested.
- `--dataset-workers=7` caps concurrent EIV–GP datasets at seven, useful when memory is limiting. The automatic cap is `floor(cores / chain-workers)`.
- `--competitor-workers=8` caps independent competitor fits at eight. Omit it or use `auto` to use the allocation.

Worker counts never exceed the supplied core allocation. Actual CPU utilization can be lower during reference calculations, prediction, disk I/O, or when fewer datasets remain in the current setting. The EIV–GP and competitor phases do not overlap. This is a process/thread budget, not OS CPU affinity. Memory use increases with concurrent datasets; the runner does not infer a safe memory budget from CPU count. Numerical-library threads are limited to one per worker by the launcher. Use terminal `Rscript` on Linux or macOS; Windows falls back to serial execution.

### Resuming an existing recovery after this upgrade

Stop the old recovery coordinator before starting the updated script, then use the **same archive, data and output paths**. Saved completed fits and scores are reused; an interrupted fit that has not reached its checkpoint must restart. The runner recognizes the immediately preceding core-budget runner by its exact source hash and permits this scheduler-only upgrade only when every other source hash, input identity, R version and dependency version is unchanged. It preserves the old provenance under `provenance-history/`. Other source changes still require a separate output folder. Changes to `--cores` or worker limits alone do not invalidate completed work.

Never run two coordinators against the same output. If the old process leaves `recovery.lock`, inspect its owner and confirm it has stopped before removing that stale lock. Worker failure is recorded without dropping another dataset's completed results; rerunning `run` resumes missing work.

The script sources the revised computation layer directly; reinstalling the package is not required to run recovery. To update a separately installed package, run `R CMD INSTALL eivGP`. If dependencies are missing, use the repository's `experiments/install_eivgp_dependencies.R` first.

## Required source layout

Under `--archive-root`, the script finds one publication archive per selected study, either directly or under `results/study1/` and `results/study2/`:

```text
study1-publication-v2-a043ddbaa2e7/
  run_summary.rds
  combined/all_raw_outputs.rds
study2-publication-v2-9ac4b52b4b1f/
  run_summary.rds
  combined/all_raw_outputs.rds
competitor-reports/
  study1/publication/{metrics,statuses,predictions}.csv
  study2/publication/{metrics,statuses,predictions}.csv
```

The ordinary Linux repository layout also works with `--archive-root=reproduction`. If multiple publication archives are present for a study, provide a root containing just the intended archives; the script refuses to guess which run to combine.

## What changes

- Oracle pool, oracle prediction, and oracle mean are separately guarded reference tasks. Failures are recorded in `reference_status`; they do not prevent EIV–GP fitting or predictive scoring. Reference computations preserve the caller's random-number state.
- The reference mean now uses two-dimensional Gauss–Hermite quadrature under the generating independent Gaussian or logistic score model. Orders 40, 80, 160, and 240 are tried until consecutive conditional means and pattern probabilities pass the declared 1e-5 absolute/relative refinement checks. Diagnostics identify quadrature explicitly; it is not presented as exact sampling or assigned a Monte Carlo SE. Failure remains visible as unavailable reference truth. The exact rejection implementation remains available through `STUDY2_ORACLE_MEAN_METHOD = "rejection"`.
- EIV–GP fits are checkpointed immediately, before predictive/mean evaluation and diagnostics. Predictive scores are checkpointed before subsequent scientific tasks. The recovery merger can retain valid predictive scores even if a later task fails. Such failures remain in `attempts.csv`; prediction completeness does not certify completion of every secondary target or convergence.
- Completed original scores are preserved. Only missing or invalid method/dataset results are retried, using the same dataset IDs. The recovered records retain separate code/source/session provenance. Historical and new reference-mean calculations may use different numerical integration methods; their method labels are retained. Table 2 does not use oracle mean truth.
- Competitor recovery uses its own cache under the output folder, avoiding legacy locks. No old locks are deleted. A locked method is reported separately without dropping other methods' results. Locks record owner PID/host/time, and cleanup only removes a lock with the matching ownership token. If an interrupted recovery leaves `recovery.lock`, inspect its owner and confirm the old process has stopped before manually removing that coordinator lock.
- UC–GP first uses the original eight starts. On failure it tries 32 and then 64 starts, with deterministic seeds offset by 1,000 per retry. Each attempt chooses the highest training likelihood among converged finite fits; recovery stops at the first valid attempt.
- EzGP first uses the original 100-evaluation budget, then 300 and 1,000 after failure. It preserves the nugget grid, training-only CV score, and CV fold seed. Invalid predictive variances are not clamped into valid ones. Every attempt is logged. The same failure-triggered rescue rule applies to every dataset; successful base fits need no rerun. LVGP retains its existing rescue protocol.

## Outputs and Table 2

The default output directory is `reproduction/recovery-publication/`:

- `provenance.rds`: source/config/code hashes and runtime information. Changing inputs, model/computation code, R/package versions, or study selection requires a new output folder; the audited scheduler-only upgrade above preserves a provenance history.
- `state.rds`: resumable combined results.
- `completeness.csv`: every method/calibration, expected IDs, missing IDs, duplicate checks, and valid counts.
- `table2_completeness.csv`: the 12 Table 2 rows.
- `predictive_summary.csv`: updated means and across-dataset MCSEs for both studies.
- `table2_complete.tex`: a `[H]` table with highlighted lowest displayed losses and no redundant R column, produced only after input checks and all 12 rows contain valid results for original IDs 1–50. Otherwise this file contains only an explicit incomplete notice.
- `study1/` and `study2/`: merged raw EIV–GP CSV/RDS outputs and competitor metrics/statuses/predictions.
- `attempts.csv`, `fits/`, and `competitor-cache/`: recovery logs and reusable checkpoints.

The supplied archive audit requires **73 Study II EIV–GP fits**, **34 Study II competitor results**, and **four Study I EzGP results**. Table 2 alone needs 25 EIV–GP fits at calibration 50 and 34 competitor results; the additional 48 EIV–GP fits complete the other calibration sizes.

The script certifies **valid predictive outputs**, not MCMC convergence. Diagnostic warnings remain included. Persistent method failures are reported, and R=50 is never forced by replacing datasets, fabricating predictions, accepting invalid variances, or changing counts. After a successful `check`, use `table2_complete.tex` to replace the current Table 2 and update the surrounding numerical discussion from `predictive_summary.csv`.

## Validation performed for this revision

Tests also verify 56-core allocation, dynamic dataset-slot refill, nested dataset/chain processes, worker-error isolation, concurrent EIV–GP result merging and resumption, restricted provenance migration, actual parallel MCMC, distinct competitor worker process IDs, identical serial/parallel predictive scores, and resumption without refitting. Tests cover forced oracle errors followed by actual tiny MCMC fitting, fit/checkpoint reuse, RNG isolation, analytic oracle checks, bounded retries, lock ownership, preservation of original successful scores, missing/duplicate ID detection, and the complete recovery orchestration with mocked expensive fits. The new reference integration passed its refinement criterion on all 25 previously failing generated datasets. On the actual cached Gaussian four-proxy UC–GP dataset 45, the original eight starts failed again and the 32-start rescue produced a valid fit. These checks do not substitute for running the remaining publication fits.
