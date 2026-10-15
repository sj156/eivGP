# Ocean input and workflow identities

- `upstream.json`: source commits, development provenance and workflow defaults.
- `inputs.csv`: row counts and SHA256 for 33 input CSVs: 16 each for the 200/500 frozen designs and one raw 1000-profile cohort.
- `validation.json`: checks executed for this integration and their scope, including historical validation records.

The installed `module-manifest.csv` records source module hashes. Output roots store package, functions and data identities, per-task settings, checksums and completion receipts.

## Current workflow (0.4.0.9007)

`ocean.R --action=fit` saves fits, analyzes saved chains using one worker, and generates HTML. `--action=continue` explicitly extends screen-selected fits and then updates analysis and HTML. Independent `analyze` and `report` actions remain available. The Rmd reads saved results; it does not fit models.

`cores` is a total fitting budget shared across fold/calibration tasks and chains. Each EIV fit receives at most its chain count; comparators receive one core. Analysis and HTML use one worker after fitting completes. Paper defaults are four chains, 100,000 iterations per chain including 5,000 burn-in; `--chains`, `--iter` and `--burn` configure a new fit. Continuation preserves initial settings and adds `--continuation-iter` iterations per existing chain (paper default 20,000). Each fit supports one explicit extension.

The main README and Ocean guide default to ocean500. The single versioned HTML illustrates a complete 200-profile five-fold short-chain run with 19 tables, 10 figures and nine Supplementary blocks. The 500-profile inputs are frozen; the 1000-profile cohort remains raw. `data/ocean_silicate/prepare_folds.R` creates fold inputs using R under `reproduction/ocean/prepared/`.

Validation covers package tests, literal main-README short-chain fit/continue commands, shared CPU reservations, serial analysis, automatic HTML, configurable initial and continuation budgets, checkpoint protection and a guarded report-only render. Scientific sampler, synthetic studies and input CSVs remain unchanged. No formal 100,000-iteration run or Linux runtime test is claimed. RAM scheduling has not been added.
