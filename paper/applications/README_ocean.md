# Ocean application

`ocean.R` is the single reproduction entrypoint. After installation, specify only
its total CPU allocation to run the default 500-profile paper experiment:

```sh
Rscript paper/applications/ocean.R --cores=16
```

With options but no `--action`, the action is `fit`: fit, analyze, export all
paper assets, then render HTML. With no arguments, print the plan. The explicit
`--action=plan|fit|analyze|continue|report` interfaces remain available. Dataset
and profile are scientific selections; worker counts are never additional knobs.

The dependency stages are:

| Stage | Tasks and dependencies | Scheduling |
|---|---|---|
| Fit | Frozen folds → 15 EIV fits, five native baseline groups, 30 CC fits | Shared pool; at most four chains per fit; independent baselines run concurrently |
| Predict | Saved fits → fixed seeded draw blocks, OI/MI scores and U intervals | Workers borrow free pool slots; fit residency is bounded separately |
| Functionals | Predictions → per-fold/calibration m/f estimates and integration checks | Outer jobs and inner draw blocks share the same pool |
| Paper assets | Scores → tables, cohort/imputation/surface figures, MCMC review | Independent figure files run concurrently |
| HTML | Completed analysis and exports → reader report | Serial assembly; no fitting |

`--cores=N` caps active computation across all these stages; it does not promise
100% utilization during I/O or serial work. Waiting coordinators release slots.
The application holds `.application.lock` through rendering to prevent concurrent writers.
The per-invocation `core-events/` logs record acquired/released slots; task logs
record completion and failures. Seeds, chains, iterations and evaluation budgets
are independent of the allocation. The CLI caps workers to detected CPU capacity.
macOS/Linux use forked workers; Windows currently falls back to serial execution.
Parallel PNG exports use Cairo instead of the fork-unsafe macOS Quartz device;
PDF uses the standard PDF device. Without Cairo on macOS, graphics export is serial.

A failed native competitor is disclosed and does not prevent analysis of valid
EIV fits. Missing/invalid EIV prerequisites and incompatible provenance still
stop execution. Analysis never retries optimizers or extends MCMC implicitly.
Rerun the same fit command to retry failed tasks and reuse compatible completed
ones; use `--action=continue` only for an intentional chain extension.

### Paper tables and figures

`reports/paper/` contains 8 named PDF figures and 6 primary LaTeX/CSV/HTML tables,
plus complete-case and paired predictive evidence. `artifact_manifest.csv` maps
every main-text/appendix label to a file, cohort, profile, availability and hash.
It includes the previously pending all-fold category-mean figure.
`method_status.csv` records missing/failed comparators independently of artifact
availability. `SUMMARY.txt` distinguishes complete output, output with method
failures, incomplete artifacts, and quick/pilot workflow checks. Completion never
certifies MCMC convergence. Raw scores and diagnostics remain beside these exports.
The notebook uses the same package computations; exports do not depend on executing
an interactive notebook. Nothing is copied into Overleaf automatically.

The current manuscript describes an older **200-profile** experiment. The default
`ocean500` regenerates the same types of tables/figures for the updated cohort;
it does not reproduce the older numerical entries. Exact historical values require
the original inputs, saved fits, selected continuations and matching environment.
Use `--dataset=ocean200` for the 200-profile design, without relabeling a new run as
the historical realization. A changed code/package identity needs a separate output
root or its original environment; never delete the saved identity to force reuse.

## Install

From the repository root, install this checkout (eivGP 0.4.0.9011):

```sh
Rscript tools/build_package.R
Rscript -e 'dir.create(".r-library",showWarnings=FALSE); install.packages("eivGP",repos=NULL,type="source",lib=".r-library")'
```

The CLI automatically finds `.r-library/`. For another installation directory,
export `EIVGP_R_LIBRARY` in the Terminal that runs the commands.
Dependencies: dplyr, tidyr, ggplot2, knitr, patchwork, posterior, TruncatedNormal,
digest, rmarkdown, kergp 0.5.8, LVGP 2.1.6 and EzGP 0.1.0. HTML requires Pandoc
on PATH or `RSTUDIO_PANDOC`; the CLI checks this before starting computation.
Restart R after updating an installed package.

## Choose the dataset

The CLI option `--dataset` and Rmd parameter `dataset` select the cohort;
both default to `ocean500`, the main application dataset.
Use `--data-root` (CLI) or `data_dir` (Rmd) for a custom prepared directory.

| Dataset | Versioned data | Train/test per fold | Training calibration | Measured test u |
|---|---|---:|---:|---:|
| `ocean200` | `data/ocean_silicate/n200` | 160/40 | 0/32/80 | 8/20 |
| `ocean500` | `data/ocean_silicate/n500` | 400/100 | 0/80/200 | 20/50 |
| `ocean1000` | `data/ocean_silicate/n1000` — raw cohort | 800/200 after preparation | 0/160/400 | 40/100 |

The 500-profile cohort is `seed_465000_joint_quintile_pressure`. The raw
1000-profile cohort is `seed_144069_oxygen_interleaved`; its fold assignments,
training/test files and masks are not versioned. See the [data guide](../../data/ocean_silicate/README.md).

## Prepare folds separately

The 200 and 500 datasets include their five-fold inputs. Prepare the 1000 dataset
once before running the experiment:

```sh
Rscript data/ocean_silicate/prepare_folds.R ocean1000
```

This writes `reproduction/ocean/prepared/ocean1000`, which both entrypoints automatically
use for `dataset="ocean1000"`. Generated files remain outside the versioned data.
The script uses base R for the complete preparation and `jsonlite` for metadata.
Grouping and balancing rules are retained. R's random numbers produce reproducible
R-specific partitions; the supplied 200/500 fold files remain frozen.

To regenerate the 500-profile design separately, or choose an output directory:

```sh
Rscript data/ocean_silicate/prepare_folds.R ocean500
Rscript data/ocean_silicate/prepare_folds.R ocean1000 /external/ocean1000-inputs
```

For custom prepared inputs, use `--data-root` in the CLI and `data_dir` in the Rmd.

**Fold method:** within each ordinal class, the 500 design forms oxygen quintiles
and pressure-ordered groups of five; the 1000 design forms successive
oxygen-ordered groups of five. Each group contributes one profile to each fold.
Assignment and within-group swaps balance standardized oxygen, log pressure,
salinity, their squares and cross-products. Training calibration is nested at
20%/50%, with class/cell coverage and boundary anchors. Test-u masks use R seed
`202609260 + fold`, with nested class-balanced 20%/50% quotas.

## Fit and explicitly continue

Run from the repository root. The main application uses ocean500.

```sh
Rscript paper/applications/ocean.R --action=plan --profile=paper --dataset=ocean500 --cores=16
Rscript paper/applications/ocean.R --action=fit --profile=paper --dataset=ocean500 --cores=16
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=20000 --cores=16
```

Results default to `<repository>/results/<dataset>/` for paper runs; pilot and quick use `results/<dataset>_pilot/` and `results/<dataset>_quick/`. Use `--output` to select any custom or existing directory. Existing results are not moved. To read an older default run, explicitly pass its `~/eivGP-ocean/<dataset>/<profile>/` path.

Sampling is fixed by profile:

| Profile | Chains | Initial iterations per chain | Burn-in |
|---|---:|---:|---:|
| paper | 4 | 50,000 | 5,000 |
| pilot | 4 | 80 | 30 |
| quick | 1 | 80 | 30 |

Initial iterations include burn-in. The CLI does not accept `--chains`, `--iter`
or `--burn`. Use `pilot` or `quick` for workflow checks, not formal results.

`cores` is a shared workflow budget for chains, prediction blocks, diagnostics,
scoring and independent figure files. Four paper chains remain fixed. Waiting
coordinators release slots, and completed chains return their slots immediately;
subsequent stages can borrow all available slots. Simultaneously resident fits
are bounded separately for memory. Analysis and report calculations retain the
requested budget; HTML assembly is serial. Windows currently uses serial
execution. `run` aliases `fit`; fitting and analysis remain separate stages.

Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1` before starting R. The pool bounds active R workers;
these startup settings bound numerical-library threads inside them. Scheduling
logs count outer job coordinators and should not be read as measured CPU usage.

Prediction uses fixed blocks and seeds independent of worker count (stream
schema 1); all prediction rows belonging to a joint draw stay together. The
new prediction stream can differ from older serial Monte Carlo realizations.
Source and cache identity checks remain strict; this update does not rewrite
existing results or migrate incompatible caches. A saved-fit `report` can use a
fresh `--diagnostics-dir` for new functional calculations without refitting.


`continue` explicitly extends initially screen-flagged fits, then updates analysis
and HTML. It preserves chain count and initial settings; burn-in is not repeated.
Each fit supports one extension, with paper default +20,000 iterations per chain
and pilot default +80. `--continuation-iter` may choose another added length:

```sh
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=10000 --cores=16
```

Passing fits are skipped; conflicting saved extension lengths are rejected.
`--folds`, `--cells`, `--modules` and `--task-ids` select tasks. `--import-dirs`
imports compatible verified tasks; `--data-root` selects prepared inputs.
New result roots use the versioned `ocean-flat-v1` layout. Scientific tasks are siblings rather than nested under method/fold/calibration wrappers:

```text
results/ocean500/
  eiv_fold01_O0/             # fit, predictions, scores and diagnostics together
  eiv_fold01_O20/
  eiv_fold01_O50/
  eiv_fold01_O20_continued/  # continuation kept distinct from initial draws
  baseline_fold01/          # fold-shared masks and original baselines
  cc_fold01_O20_EzGP/        # one native competitor task
  data/                     # frozen input snapshot
  logs/paper_fit/            # execution status and task logs
  reports/diagnostics/      # derived tables, figures and caches
  reports/prediction_report/
  _meta/                    # original provenance, execution settings
    storage-layout.rds      # physical layout schema, not scientific identity
    tasks/                  # task configuration and completion receipts
  ocean.html
```

Other folds follow the same naming. Pilot/quick roots remain separate from paper roots. `_meta/initial` and `_meta/continuation` retain any stage-wide legacy metadata; task payloads never live beneath them. Explicit report directories remain supported.

Unmarked existing roots retain their legacy paths. Readers, reports, continuation, completion receipts and imports resolve through the package's storage-path API. Importing between layouts translates paths and checks payload hashes. Scientific provenance checks are unchanged: migration does **not** make a historical fit compatible with a different sampler/package version. Derived diagnostic caches still require matching code and settings; rebuild them in a new diagnostics directory if those changed.

### Plan migration after the existing job finishes

The migration script is **dry-run by default** and refuses a source with `.workflow.lock`. Do not remove a lock while its job is running. Use a new, absent destination; existing outputs are never overwritten.

```bash
Rscript paper/applications/migrate_ocean.R \
  --source=/path/to/results/ocean/ocean500/paper \
  --destination=/path/to/results/ocean500
```

This prints a collision-checked file map without copying or moving results. After reviewing the plan and confirming the source job has finished, explicit `--apply=true` copies into a temporary destination, verifies source/destination hashes, rewrites receipt paths only, and publishes the completed directory. The original source is retained. It needs enough free space for a second copy. `_meta/migration-manifest.csv` records both sets of hashes; `_meta/migration.rds` retains the source identity. Existing rendered HTML and provenance text remain historical artifacts; regenerate reports through the updated reader when needed. No scientific draws are thinned, deleted, recomputed or relabeled.

When preparing this update alongside an active job, install into a separate R library and select it with `EIVGP_R_LIBRARY`; do not replace the library used by the active process.

## Generate HTML from saved results

`fit` and `continue` already generate HTML. To update results without MCMC:

```sh
Rscript paper/applications/ocean.R --action=analyze --profile=paper --dataset=ocean500
Rscript paper/applications/ocean.R --action=report --profile=paper --dataset=ocean500
```

`analyze` reads saved chains and rebuilds predictions/diagnostics/HTML; `report`
reads saved predictions and renders HTML. Direct Rmd rendering is also available:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean500",profile="paper",action="report",
    output_dir="",
    run_id="ocean500_paper_report",appendix=TRUE))
```

Use `analyze` to create predictions from any checkpoints produced through the low-level fitting API.
The Rmd defaults to `report`, also accepts `plan`, and rejects fitting actions.
Custom package libraries must also be on the R session's `.libPaths()`.

`paper/config/ocean200.yml`, `ocean500.yml` and `ocean1000.yml` are report parameter
templates. Read one with `yaml::read_yaml()` and pass it as `params`.

OI/CC/MI protocols and scoring scopes are defined in the report. Masks, counts
and eligible best-score emphasis are computed from inputs/results. Keep fitted
objects and chains in external output roots. Input hashes and validation scope
are recorded in `paper/manifests/ocean/`.
