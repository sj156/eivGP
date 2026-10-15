# Ocean application

`ocean.R` is the command-line experiment entrypoint: it fits/resumes models and
saves predictions, summaries and diagnostics. `ocean.Rmd` reads saved results and
generates the reader-facing HTML; it does not fit models or continue MCMC. Only `ocean.html` is versioned as the complete five-fold, short-chain
`pilot` example report for the 200-profile cohort.
The package follows the upstream structure: Ocean execution and core model calls
are in `helpers.R`; tables, figures and diagnostics are in `plotting.R`.
The shared `mcmc.R` sampler is unchanged.

## Install

From the repository root, install this checkout (eivGP 0.4.0.9003):

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

## Run experiments from the command line

Run from the repository root after installation. The runner also finds the
checkout from its own path when launched elsewhere. Set `EIVGP_R_LIBRARY` when
using a custom package library; otherwise `.r-library/` is used when present.
The installed source identity is checked against this checkout.

```sh
# If installed in library/ above, set this in the Terminal session:
export EIVGP_R_LIBRARY="$(pwd)/library"

# Validate inputs and list tasks; no fitting. This is also the no-argument default.
Rscript paper/applications/ocean.R --action=plan --profile=paper --dataset=ocean500

# Minimal fitting check: fold 1, one chain, 80 iterations / 30 burn-in.
Rscript paper/applications/ocean.R --action=run --profile=quick --dataset=ocean500 --cores=1

# All five folds and report modules, using short chains.
Rscript paper/applications/ocean.R --action=run --profile=pilot --dataset=ocean500 --cores=4 --diagnostics-workers=2

# Complete experiment with formal numerical settings.
Rscript paper/applications/ocean.R --action=run --profile=paper --dataset=ocean500 --cores=4 --diagnostics-workers=2
```

Without `--output`, results are written to `results/ocean/<dataset>/<profile>/`.
Use `--output=/external/ocean500-paper` to choose another root. `--data-root`
overrides the prepared input directory. Run `--help` for task selectors and
separate diagnostics/import directories. Existing parameter RDS files remain
accepted as the runner's single positional argument.

`quick` covers EIV/CC prediction only; it omits OI comparators, MCMC review,
response-function plots and input-imputation plots. `pilot` covers all modules:
four chains, 80 initial iterations, 30 burn-in and one diagnostic-selected
80-iteration continuation. `paper` uses four chains, 20,000 initial iterations,
5,000 burn-in and one diagnostic-selected 20,000-iteration continuation.
Short-chain completion validates execution, not posterior convergence.

Reissue the same command to resume compatible completed tasks. Changes to inputs,
package source or numerical settings require a fresh output root. Select
`--folds=1,3` or `--folds=1:5`, `--cells=O20,O50`, `--modules=...` or
`--task-ids=...` to distribute work; keep separate output roots per process.
`--import-dirs=/path/one,/path/two` combines completed compatible tasks.

## Generate HTML from saved results

After the runner finishes, render the Rmd in R with the same dataset, profile and
output root. No fitting or continuation is performed by the notebook.

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean500",profile="paper",action="report",
    output_dir="results/ocean/ocean500/paper",diagnostics_workers=2,
    run_id="ocean500_paper_report",appendix=TRUE))
```

For short-chain results, use the same run profile (`quick` or `pilot`) and its
corresponding output root. When the CLI used a custom `--output`, give that same
directory to the Rmd's `output_dir`. If task selection
was restricted, retain the same `folds`, `cells` and `modules` in the report.
An installed package in a custom library must also be on the R session's
`.libPaths()`. HTML rendering requires Pandoc.

The Rmd defaults to `action="report"`. `action="plan"` can display the proposed
analysis without fitting; `action="run"` is rejected with the CLI command to use.
`ocean.R --action=report` rebuilds saved-result summaries/figures without fitting;
it does not render HTML. Use the Rmd to produce HTML.

`paper/config/ocean200.yml`, `ocean500.yml` and `ocean1000.yml` are Rmd report
parameter templates. Read one with `yaml::read_yaml()` and pass it as `params`.

OI/CC/MI protocols and scoring scopes are defined in the report. Masks, counts
and eligible best-score emphasis are computed from inputs/results. Keep fitted
objects and chains in external output roots; use a fresh root when inputs,
package identity or numerical settings change. Input hashes and validation scope
are recorded in `paper/manifests/ocean/`.
