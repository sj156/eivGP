# Ocean application

`ocean.R` runs fitting, single-worker analysis and HTML generation in one command.
Explicit continuation updates the saved chains, analysis and HTML. Internally,
fitting and analysis remain separate stages. `ocean.Rmd` is the reader-facing
report and reads saved results; it cannot fit models or continue MCMC.
Only `paper/applications/ocean.html` is versioned as a complete example report.
Ocean execution lives in `helpers.R`, reporting in `plotting.R`; the shared
`mcmc.R` sampler is unchanged.

## Install

From the repository root, install this checkout (eivGP 0.4.0.9007):

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

Run from the repository root. These commands default to the 500-profile cohort.

```sh
# Inspect the tasks and resource plan without fitting.
Rscript paper/applications/ocean.R --action=plan --profile=paper --dataset=ocean500 --cores=16 --chains=4 --iter=100000

# Fit, then analyze using one worker and generate HTML.
Rscript paper/applications/ocean.R --action=fit --profile=paper --dataset=ocean500 --cores=16 --chains=4 --iter=100000

# Extend initially flagged fits, then update analysis and HTML automatically.
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=20000 --cores=16
```

`--cores` is a total fitting budget. A shared queue runs independent
fold/calibration tasks and assigns up to `--chains` cores to each EIV fit;
comparator tasks receive one core each. All active allocations together stay
within the total budget. Completed tasks free resources for the next task.
The four-chain examples can use 12, 16 or 56 cores for 3, 4 or 14 simultaneous
EIV fits when enough work is pending. More chains can be selected explicitly
before fitting; the resource planner never changes the statistical chain count.
Windows falls back to serial execution.

After the fitting queue finishes, the coordinator analyzes saved chains using
one worker and renders `ocean.html` in the result root. There is no analysis-worker
option or overlap between fitting and analysis. Low-level package fit/analyze
stages remain separately callable. `run` aliases `fit`.

Defaults: paper uses four chains, 100,000 initial iterations per chain and 5,000
burn-in; pilot uses four chains, 80 initial iterations and 30 burn-in. `quick`
uses one chain and fold 1. Set `--chains` and `--iter` before starting a new
paper/pilot fit; `--burn` is optional. `iter` counts total iterations per chain,
including burn-in, and must exceed `burn`. For example, a new fit can use:

```sh
Rscript paper/applications/ocean.R --action=fit --profile=paper --dataset=ocean500 --chains=8 --iter=60000 --burn=5000 --cores=16 --output=results/ocean/ocean500/custom
```

Later commands infer the saved initial chain count, iteration count and burn-in.
`continue` retains these settings: requesting a different `--chains`, `--iter`
or `--burn` is rejected. Use `--continuation-iter` for **additional iterations
per existing chain**; burn-in is not repeated. `cores` can change between stages.

`continue` performs at most one explicit extension on fits flagged by the initial
finite-Rhat/bulk-ESS screen; it never runs automatically. Passing fits are skipped.
The default added iterations per chain are paper 20,000 and pilot 80. Choose a
different length at the continuation step:

```sh
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=10000 --cores=16
```

Saved extension lengths drive reuse and report labels. Requesting a different
length for an already-started extension fails. Verified checkpoints survive
analysis/report failures; rerun `analyze` or `report` after correcting the issue.

Results default to `results/ocean/<dataset>/<profile>/`; `--output` overrides it.
`--folds`, `--cells`, `--modules` and `--task-ids` select tasks. Keep one owner per
output root; changed source, data or initial-fit numerical settings require a
fresh root. `--import-dirs` imports compatible verified tasks; `--data-root`
selects custom prepared inputs. Task logs are under `jobs/<run_id>/tasks/`;
`task_status.csv` and `execution_schedule.csv` retain PIDs, core allocations,
completion and failures.

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
    output_dir="results/ocean/ocean500/paper",
    run_id="ocean500_paper_report",appendix=TRUE))
```

Use `analyze` to create predictions from any checkpoints produced through the low-level fitting API.
The Rmd defaults to `report`, also accepts `plan`, and rejects fitting actions.
Custom package libraries must also be on the R session's `.libPaths()`.

`paper/config/ocean200.yml`, `ocean500.yml` and `ocean1000.yml` are report parameter
templates. Read one with `yaml::read_yaml()` and pass it as `params`.

OI/CC/MI protocols and scoring scopes are defined in the report. Masks, counts
and eligible best-score emphasis are computed from inputs/results. Keep fitted
objects and chains in external output roots; use a fresh root when inputs,
package identity or numerical settings change. Input hashes and validation scope
are recorded in `paper/manifests/ocean/`.
