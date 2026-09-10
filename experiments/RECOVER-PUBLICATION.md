# Recover the original 50-dataset publication comparisons

Run this from the repository containing the revised `codes/`, `eivGP/`, and `experiments/` directories. Use the Linux machine with the **original frozen `.rds` inputs**. Do not regenerate datasets to replace failures. Recovery does not modify the archives, original reports, or the Overleaf paper.

## Commands

From the repository root on Linux, with the usual `reproduction/` layout:

```sh
Rscript --vanilla experiments/recover_publication.R plan
Rscript --vanilla experiments/recover_publication.R run > recovery.log 2>&1
Rscript --vanilla experiments/recover_publication.R check
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

`--chain-workers=4` is the default: one dataset runs at a time, with up to four chains in parallel. This preserves the archived iterations, warmup, priors, seeds, and calibration subsets. The code does not estimate runtime or reduce the MCMC budget. Keep the output directory and its fit checkpoints for resumption.

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

- `provenance.rds`: source/config/code hashes and runtime information. Changing inputs, code, R/package versions, or study selection requires a new output folder.
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

Tests cover forced oracle errors followed by actual tiny MCMC fitting, fit/checkpoint reuse, RNG isolation, analytic oracle checks, bounded retries, lock ownership, preservation of original successful scores, missing/duplicate ID detection, and the complete recovery orchestration with mocked expensive fits. The new reference integration passed its refinement criterion on all 25 previously failing generated datasets. On the actual cached Gaussian four-proxy UC–GP dataset 45, the original eight starts failed again and the 32-start rescue produced a valid fit. These checks do not substitute for running the remaining publication fits.
