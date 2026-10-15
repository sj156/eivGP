# Ocean application

`ocean.R` is the command-line experiment entrypoint: it separates fitting, explicit continuation and single-worker analysis. `ocean.Rmd` reads saved results and
generates the reader-facing HTML; it does not fit models or continue MCMC. Only `ocean.html` is versioned as the complete five-fold, short-chain
`pilot` example report for the 200-profile cohort.
The package follows the upstream structure: Ocean execution and core model calls
are in `helpers.R`; tables, figures and diagnostics are in `plotting.R`.
The shared `mcmc.R` sampler is unchanged.

## Install

From the repository root, install this checkout (eivGP 0.4.0.9004):

```r
dir.create("library",showWarnings=FALSE)
.libPaths(c(normalizePath("library"),.libPaths()))
source("tools/build_package.R")
install.packages("eivGP",repos=NULL,type="source",lib=.libPaths()[1])
# Let command-line runners find this same installation:
Sys.setenv(EIVGP_R_LIBRARY=normalizePath("library"))
```

Dependencies: dplyr, tidyr, ggplot2, knitr, patchwork, posterior, TruncatedNormal,
digest, rmarkdown, kergp 0.5.8, LVGP 2.1.6 and EzGP 0.1.0. HTML rendering requires
Pandoc. Restart R after updating an installed package.

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

## Fit, analyze and explicitly continue

Run from the repository root after installation. Set `EIVGP_R_LIBRARY` to a
custom installation directory; otherwise the CLI uses `.r-library/` when present.
The installed module identity must match this checkout.

```sh
# If installed in library/ using the instructions above:
export EIVGP_R_LIBRARY="$(pwd)/library"

# Validate inputs and inspect the task catalogue.
Rscript paper/applications/ocean.R --action=plan --profile=paper --dataset=ocean500

# Fit initial models and save checkpoints; no automatic continuation or EIV analysis.
Rscript paper/applications/ocean.R --action=fit --profile=paper --dataset=ocean500 --cores=4

# Analyze saved chains, predict, score and plot; always one worker, no MCMC.
Rscript paper/applications/ocean.R --action=analyze --profile=paper --dataset=ocean500

# If the initial diagnostic screen flags fits, explicitly request one extension.
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --cores=4

# Analyze the newly extended chains before rendering the final report.
Rscript paper/applications/ocean.R --action=analyze --profile=paper --dataset=ocean500
```

`fit` saves EIV checkpoints and package-generated basic diagnostics. Native
comparators are fitted through their existing package adapters, which also
return native prediction outputs. EIV predictions, imputations, scores and full
report diagnostics belong to `analyze`. `run` is an alias for `fit`.

`continue` performs at most the existing one-round extension on selected fits
flagged by the initial finite-Rhat/bulk-ESS screen; it never runs automatically.
A fit passing that screen is skipped. Repeating the same stage reuses compatible
outputs. Initial results remain available for analysis before an extension is
requested; convergence flags stay visible.

The `paper` profile retains four chains, 20,000 initial iterations, 5,000 burn-in
and a 20,000-iteration explicit extension. `pilot` retains four chains, 80 initial
iterations, 30 burn-in and an 80-iteration extension. `quick` uses one fold and
one chain; its task catalogue omits continuation and the full diagnostic modules.
Analysis is always serial; there is no diagnostics-worker option. `cores` controls
MCMC chain parallelism, not the number of chains.

Results default to `results/ocean/<dataset>/<profile>/`; `--output` overrides it.
Keep the same dataset, profile, output root and task selections across stages.
Use a fresh root for changed inputs, installed source or numerical settings.
`--folds=1,3`, `--cells=O20,O50`, `--modules=...` and `--task-ids=...` select work.
Give each process its own output root. `--import-dirs` imports only compatible
completed tasks; `--data-root` selects custom prepared inputs. See `--help`.

## Generate HTML from saved results

After `analyze`, render the reader-facing Rmd using the same inputs and result
root. It reads saved fits/results and cached calculations; it cannot fit or extend
MCMC. Any necessary report-side calculation also runs in one worker.

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean500",profile="paper",action="report",
    output_dir="results/ocean/ocean500/paper",
    run_id="ocean500_paper_report",appendix=TRUE))
```

`ocean.R --action=report` reads existing predictions and rebuilds report summaries
without fitting; use `analyze` first to create predictions from new checkpoints.
The Rmd defaults to `report`, also accepts `plan`, and rejects fitting actions.
Custom package libraries must also be on the R session's `.libPaths()`.

`paper/config/ocean200.yml`, `ocean500.yml` and `ocean1000.yml` are report parameter
templates. Read one with `yaml::read_yaml()` and pass it as `params`.

OI/CC/MI protocols and scoring scopes are defined in the report. Masks, counts
and eligible best-score emphasis are computed from inputs/results. Keep fitted
objects and chains in external output roots; use a fresh root when inputs,
package identity or numerical settings change. Input hashes and validation scope
are recorded in `paper/manifests/ocean/`.
