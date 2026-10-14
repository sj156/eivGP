# Independent Ocean application layer

This application calls the unmodified upstream eivGP 0.4.0 package at commit `cac1498d6794ece85fbcb2b47d17ad26fdce4ddd`. Its three package-module hashes are pinned in `UPSTREAM_PACKAGE.csv`. Source changes are confined to this application layer; it does not call the upstream `run_ocean_analysis()` workflow.

## Layout

- `paper/applications/ocean.R`: shared Rmd/batch entry; checks the pinned package, executes selected tasks and assembles available results.
- `paper/applications/ocean.Rmd`: scientific setup, interpretation, and OI/CC/MI result displays.
- `paper/applications/ocean.html`: generated report. Inspect its computation profile/fold counts before treating it as a paper result.
- `paper/applications/ocean_helpers/01_setup.R`: masks, selectors, receipts and execution.
- `02_fit_models.R`: application orchestration; calls upstream fit/predict/impute/continue and native competitor adapters.
- `03_reporting.R`: result reading, labels, paired/pool comparisons, automatic ranking and tables.
- `04_diagnostics.R`: numerical review, imputation, mean/surface plots and empirical-reference histograms.
- `data/ocean_silicate/`: the frozen cohort, ten train/test CSVs and five explicit test-mask CSVs.

The discarded observation-status/zero-fill comparator has no application selector, encoder, fit function or result reader here. The upstream package is left intact. Old historical fit archives are retained separately and are never relabeled as CC.

## Package environment

Use a fresh R session and a project-specific package library. Install the fixed upstream commit, with required dependencies prepared first:

```r
dir.create("library",showWarnings=FALSE)
.libPaths(c(normalizePath("library"),.libPaths()))
remotes::install_github("sj156/eivGP",subdir="eivGP",
  ref="cac1498d6794ece85fbcb2b47d17ad26fdce4ddd",
  lib=.libPaths()[1],upgrade="never")
```

Required packages include kergp 0.5.8, LVGP 2.1.6, EzGP 0.1.0, digest, posterior, TruncatedNormal, dplyr, tidyr, ggplot2, patchwork, knitr and rmarkdown. Pandoc is required for HTML. Neither sourcing nor rendering installs dependencies. Do not overwrite an existing 0.3.1 library or load two eivGP versions in one R session.

## Run and report

Run from the repository root. For this independent Ocean extension, use the commands here. The old positional `Rscript ocean.R DATA_DIR OUTPUT_DIR ACTION` interface is not used by this entrypoint. The default document action is `plan`:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  output_file="ocean.html",params=list(action="plan"))
```

For a fresh execution of the entire Rmd workflow with shorter numerical budgets:

```r
rmarkdown::render("paper/applications/ocean.Rmd",output_file="ocean.html",
  params=list(profile="pilot",action="run",
    output_dir="reproduction/ocean/full_workflow_pilot_v040",
    cores=4,diagnostics_workers=2,run_id="full_workflow"))
```

The `pilot` profile retains all five folds and all three calibration levels: 15 initial EIV fits, five OI comparator groups (three methods per group), 30 CC fits, and up to 15 diagnostic-selected continuations. It also computes prediction tables, imputation diagnostics/distributions, m/f plots and saved-chain checks. Its budgets are 80 iterations, 30 burn-in, four chains, and one additional 80-iteration continuation when selected. Native comparator controls use fewer optimizer starts/iterations; function plots use eight posterior states per chain and 16 latent integration nodes. These settings verify execution, not scientific precision. Failed fitting tasks stop rendering; a completed task can still have flagged optimizer/MCMC diagnostics, which remain in the report.

For the smaller, single-fold EIV/CC-only execution:

```r
rmarkdown::render("paper/applications/ocean.Rmd",output_file="ocean.html",
  params=list(profile="quick",action="run",folds=1,
    cells=c("O0","O20","O50"),modules=c("eiv","complete_case_baselines","prediction"),
    output_dir="reproduction/ocean/quick_v040",cores=1,run_id="quick"))
```

This selection exercises three EIV conditions and six native CC fits; it does not exercise the original full-training OI baseline group or the full five-fold diagnostics/surface stage. The `pilot` and `paper` profiles include those groups. Formal runs use `profile="paper"`, a fresh result root and the existing fixed scientific settings, with collaborators selecting folds/modules/task IDs as appropriate.

For the collaborator's full-budget execution from a fresh result root:

```r
rmarkdown::render("paper/applications/ocean.Rmd",output_file="ocean.html",
  params=list(profile="paper",action="run",
    output_dir="/external/new-ocean-paper-run",
    cores=4,diagnostics_workers=2,run_id="paper_full"))
```

Replace `/external/new-ocean-paper-run` with a writable location on the collaborator's machine and choose the CPU budget explicitly. Keep the default five folds, three calibration conditions and all modules for the complete report. The default Rmd action is `plan`, so omitting `action="run"` does not fit models. Full EIV settings are four chains, 20,000 initial iterations, 5,000 burn-in, and one diagnostic-selected 20,000-iteration continuation; completing the budget does not imply that every diagnostic screen passes.

For saved results, choose `action="report"` with the actual EIV/OI result root; optionally supply `cc_output_dir` for independent CC results. Report does not refit:

```r
rmarkdown::render("paper/applications/ocean.Rmd",output_file="ocean.html",
  params=list(profile="paper",action="report",
    output_dir="/external/eiv_oi_results",cc_output_dir="/external/cc_results",
    modules=c("prediction","diagnostics","surface_comparison"),run_id="report"))
```

For batch execution, save parameter fields to an RDS and run `Rscript paper/applications/ocean.R parameters.rds`. Supported fields include profile/action/folds/cells/modules/task_ids/cores/run_id, data_dir/output_dir, import_dirs and cc_output_dir. The shared entry returns computed result objects; the Rmd supplies their scientific presentation. Source the entry and load helpers in an R session for custom dispatch.

Keep large fits, chains, reports' reconstruction caches and logs in external result directories. Package/version/application identity changes require new fitting roots; historical receipts must not be rewritten to force reuse. OI/CC/MI have the meanings stated in the report; CC training uses only measured-u profiles, and its comparison uses the identical observed-u test subset for all methods.