# Ocean application

`ocean.R` is the batch entrypoint; `ocean.Rmd` runs the analysis and presents its
results. Only `ocean.html` is versioned as the 200-profile, single-fold `quick` example report.
The package follows the upstream structure: Ocean execution and core model calls
are in `helpers.R`; tables, figures and diagnostics are in `plotting.R`.
The shared `mcmc.R` sampler is unchanged.

## Install

From the repository root, install this checkout (eivGP 0.4.0.9002):

```r
dir.create("library",showWarnings=FALSE)
.libPaths(c(normalizePath("library"),.libPaths()))
source("tools/build_package.R")
install.packages("eivGP",repos=NULL,type="source",lib=.libPaths()[1])
```

Dependencies: dplyr, tidyr, ggplot2, knitr, patchwork, posterior, TruncatedNormal,
digest, rmarkdown, kergp 0.5.8, LVGP 2.1.6 and EzGP 0.1.0. HTML rendering requires
Pandoc. Restart R after updating an installed package.

## Choose the dataset

The Rmd parameter `dataset` selects the cohort; `data_dir` can specify a prepared
directory explicitly.

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
once before running the Rmd:

```sh
Rscript data/ocean_silicate/prepare_folds.R ocean1000
```

This writes `reproduction/ocean/prepared/ocean1000`, which the Rmd automatically
uses for `dataset="ocean1000"`. Generated files remain outside the versioned data.
The script uses base R for the complete preparation and `jsonlite` for metadata.
Grouping and balancing rules are retained. R's random numbers produce reproducible
R-specific partitions; the supplied 200/500 fold files remain frozen.

To regenerate the 500-profile design separately, or choose an output directory:

```sh
Rscript data/ocean_silicate/prepare_folds.R ocean500
Rscript data/ocean_silicate/prepare_folds.R ocean1000 /external/ocean1000-inputs
```

For a custom output directory, pass it as the Rmd's `data_dir`.

**Fold method:** within each ordinal class, the 500 design forms oxygen quintiles
and pressure-ordered groups of five; the 1000 design forms successive
oxygen-ordered groups of five. Each group contributes one profile to each fold.
Assignment and within-group swaps balance standardized oxygen, log pressure,
salinity, their squares and cross-products. Training calibration is nested at
20%/50%, with class/cell coverage and boundary anchors. Test-u masks use R seed
`202609260 + fold`, with nested class-balanced 20%/50% quotas.

## Plan, run and report

The default action is `plan`, which validates inputs and lists tasks:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean200",action="plan"))
```

A minimal fitting check uses the 200 dataset and the `quick` profile:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean200",profile="quick",action="run",
    output_dir="/external/ocean200-quick",cores=1,run_id="ocean200_quick"))
```

`quick` uses fold 1, one chain, 80 iterations and 30 burn-in, covering the EIV/CC
prediction workflow. `pilot` uses all five folds and all modules with four short
chains. Full-budget fitting uses `paper`:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean1000",profile="paper",action="run",
    output_dir="/external/ocean1000-paper",cores=4,diagnostics_workers=2,
    run_id="ocean1000_paper"))
```

Formal settings retain four chains, 20,000 initial iterations, 5,000 burn-in and
one diagnostic-selected 20,000-iteration continuation. To read saved results
without refitting, use the original run profile and output directory:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean1000",profile="paper",action="report",
    output_dir="/external/ocean1000-paper",
    diagnostics_dir="/external/ocean1000-report",run_id="report"))
```

`paper/config/ocean200.yml`, `ocean500.yml` and `ocean1000.yml` provide parameter
templates. Read one with `yaml::read_yaml()` and pass it as `params`. The batch
`ocean.R` also accepts an RDS containing these fields. Select folds/cells/modules
to distribute work; use `import_dirs` to combine completed compatible tasks.

OI/CC/MI protocols and scoring scopes are defined in the report. Masks, counts
and eligible best-score emphasis are computed from inputs/results. Keep fitted
objects and chains in external output roots; use a fresh root when inputs,
package identity or numerical settings change. Input hashes and validation scope
are recorded in `paper/manifests/ocean/`.
