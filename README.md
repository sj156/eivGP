# eivGP

`eivGP` implements errors-in-variables Gaussian process regression with
numeric predictors, ordinal proxies for latent continuous inputs, and sparse
calibration measurements. This repository is both the literate package source
and the reproducibility companion for the numerical studies.

## Install

Install the generated package subdirectory from GitHub:

```r
pak::pak("sj156/eivGP/eivGP", dependencies = TRUE)
```

From a checkout, regenerate and check the package with:

```r
litr::render("index.Rmd")
devtools::check("eivGP", document = FALSE)
```

## Stable package interface

```r
fit <- eivGP::fit_eivgp(
  X, y, C, U_obs,
  engine = "multivariate",
  latent_dim = 2L,
  kernel = "matern"
)

mean_draws <- eivGP::predict_eivgp(
  fit, new_X = X_new, new_C = C_new, target = "mean"
)
response_draws <- predict(
  fit, new_X = X_new, new_C = C_new, target = "response"
)
latent_draws <- eivGP::impute_eivgp(fit, new_C = C_new)
```

The published mixed-input competitors are available through
`fit_ucgp()`, `fit_lvgp()`, `fit_ezgp()`, and
`fit_mixedgp_competitor()` when their suggested packages are installed.

## Reproduce the paper's numerical experiments

This section is the supported route for reproducing the numerical tables and
figures in the paper. Real-data analyses are intentionally not included yet.

### 1. Obtain the reproducibility companion and install its dependencies

Clone this repository, then start R from the repository root. The numerical
studies require `rmarkdown`, `pkgload`, and the competitor packages identified
by the publication preflight. Install the package dependencies with:

```r
pak::pak(c(
  "rmarkdown", "pkgload",
  "sj156/eivGP/eivGP"
))
```

Then install the package from that same checkout (repeat this after pulling
changes):

```sh
R CMD INSTALL eivGP
```

The documents use `library(eivGP)`, as an external user would. The package
installs the frozen simulation source used by its public functions. If a
required competitor is unavailable, the preflight reports the exact missing
package; a publication run stops rather than silently changing the paper's
comparison set.

### 2. Generate frozen synthetic data once

Generate each study's deterministic input data before fitting. Open
[`experiments/study1_synthetic_data.Rmd`](experiments/study1_synthetic_data.Rmd)
or [`experiments/study2_synthetic_data.Rmd`](experiments/study2_synthetic_data.Rmd),
set `run_mode: "publication"`, and render it. Each creates the revised frozen
design below the chosen reproduction root and refuses to overwrite existing files.

The same data-only operation is available from the command line:

```sh
EIVGP_ARTIFACT_ROOT=/path/to/output-directory \
EIVGP_RUN_MODE=publication \
Rscript experiments/run_publication_study.R study1-data
```

Use `study2-data` for Study II. The generated manifests and checksums are the
fixed inputs for every later fitting run.

### 3. Common manual publication procedure

Use this sequence for any ready study after its frozen inputs have been
generated. Run every command from the repository root. Do not regenerate data
when reproducing an existing study.

```sh
R CMD INSTALL eivGP
```

Select the study-specific input and output directories:

```sh
STUDY=study1
DATA_ROOT="reproduction/data/synthetic/$STUDY"
RESULT_ROOT="reproduction/results/$STUDY"
export STUDY DATA_ROOT RESULT_ROOT
```

For Study II, change the first line to `STUDY=study2`. Confirm that the frozen
cell manifests are present:

```sh
find "$DATA_ROOT" -name manifest.rds
```

Run the read-only preflight. It must show every design cell with
`complete_count TRUE` and all required packages as available:

```sh
Rscript --vanilla -e 'library(eivGP); study <- Sys.getenv("STUDY"); data_root <- Sys.getenv("DATA_ROOT"); result_root <- Sys.getenv("RESULT_ROOT"); config <- if (identical(study, "study1")) study1_simulation_config(mode = "dry_run", workers = 8L, data_root = data_root, output_root = result_root) else study2_simulation_config(mode = "dry_run", workers = 8L, data_root = data_root, output_root = result_root); if (identical(study, "study1")) run_study1_simulation(config) else run_study2_simulation(config)'
```

Start the formal fit and aggregation:

```sh
EIVGP_ARTIFACT_ROOT="$PWD/reproduction" \
EIVGP_RUN_MODE=publication \
EIVGP_WORKERS=8 \
Rscript --vanilla experiments/run_publication_study.R "$STUDY"
```

Keep the terminal open until the command exits successfully. If the machine
does not have enough memory for eight workers, restart the final command with
`EIVGP_WORKERS=4`. Do not treat smoke-mode output as publication output.

The real-data application will follow this same pattern—freeze an audited data
bundle, run a preflight, fit through `eivGP`, then write results and a report.
Its data bundle and run wrapper are not released yet, so
[`applications/real_data_application.Rmd`](applications/real_data_application.Rmd)
is intentionally not a publication command at this stage.

### Monitor a long-running fit

The runner writes a single, atomically updated progress record at
`config/progress.csv` inside its fingerprinted run directory. From a second
terminal, locate the newest Study I run and print its current phase:

```sh
RUN_DIR="$(find reproduction/results/study1 -maxdepth 1 -type d -name 'study1-publication-*' | sort | tail -n 1)"
cat "$RUN_DIR/config/progress.csv"
```

For a live terminal monitor that refreshes every 30 seconds:

```sh
while true; do clear; date; cat "$RUN_DIR/config/progress.csv"; sleep 30; done
```

`config/cell_status.csv` records completed cells and elapsed times;
`config/diagnostic_gates.csv` records their quality gates. Stop the monitor
with `Ctrl-C`; this does not stop the experiment.

### 4. Run another numerical study

Use the portable launcher from the repository root. Substitute any writable
directory for `/path/to/output-directory`; it is the only location that will
receive generated data, simulation results, tables, figures, and HTML reports.

```sh
EIVGP_ARTIFACT_ROOT=/path/to/output-directory \
  Rscript experiments/run_publication_study.R study1
```

Use `study2` in place of `study1` to reproduce Study II. The command defaults
to the frozen, reportable `publication` design. It verifies frozen inputs and
then performs these stages in order:

1. Runs the prespecified `eivGP` and competitor analyses.
2. Validates the paired design and the exact competitor implementations,
   aggregates results, and writes the paper-ready result artifacts.
3. Renders an HTML audit report.

Every run is fingerprinted. Re-run the same command after an interruption;
the frozen inputs are checksum-verified before work begins.

### 5. Locate and inspect the results

For Study I, the command writes these files beneath the chosen reproduction
directory:

```
data/synthetic/study1/                  # frozen synthetic datasets and manifest
results/study1/<fingerprinted-run>/
  cells/<design-cell>/figures/           # paper figures for each design cell
  cells/<design-cell>/tables/            # CSV and LaTeX tables for each design cell
  combined/                             # aggregated raw results and paired comparisons
  config/                               # validation records, hashes, versions, session details
reports/study1_numerical_experiment.html
```

Study II uses the analogous `study2` directories. The HTML report lists the
exact aggregated files produced for that run. Keep the full fingerprinted run
directory with any manuscript result: it is the audit trail linking tables and
figures to input data, code hashes, package versions, diagnostics, and
validation outcomes.

### Dry runs and computing controls

Before committing compute time, verify the local setup without modifying
outputs:

```sh
EIVGP_RUN_MODE=dry_run \
Rscript experiments/run_publication_study.R study1
```

For a small, explicitly non-reportable end-to-end test, use
`EIVGP_RUN_MODE=smoke`. Generate its separate smoke inputs before fitting:

```sh
EIVGP_ARTIFACT_ROOT=/path/to/output-directory EIVGP_RUN_MODE=smoke \
  Rscript experiments/run_publication_study.R study1-data
EIVGP_ARTIFACT_ROOT=/path/to/output-directory EIVGP_RUN_MODE=smoke \
  Rscript experiments/run_publication_study.R study1
```

Smoke data are stored under `data/smoke/` and are never paper results. Other
optional controls are:

- `EIVGP_WORKERS=8` sets the number of replication workers.
- `MIXEDGP_R_LIBRARY=/path/to/R-library` adds a custom R library before the
  preflight, useful on shared systems.
- `EIVGP_ARTIFACT_ROOT` is optional; without it, reproducibility files go to
  `reproduction/` under the repository checkout.

No personal filesystem path is embedded in these commands. The standalone
documents—[`experiments/study1_numerical_experiment.Rmd`](experiments/study1_numerical_experiment.Rmd)
and [`experiments/study2_numerical_experiment.Rmd`](experiments/study2_numerical_experiment.Rmd)—may also be rendered from RStudio after setting their `run_mode` parameter to `publication`.

### Real-data application

The real-data workflow will be released separately after its inputs and
analysis specification have been audited. Its placeholder document,
[`applications/real_data_application.Rmd`](applications/real_data_application.Rmd),
is not part of the current paper-reproduction instructions.

## Source organization

The cleaned R scripts under `script/R-scripts/` are the computational source
of truth. The explicitly ordered `.Rmd` chapters form the `litr` package
skeleton and copy those scripts into the generated `eivGP/` package. Historical
chapters are retained under `legacy/` for provenance and are excluded from the
active package build.
