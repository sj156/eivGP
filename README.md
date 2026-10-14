# eivGP

Ocean workflows and the 200/500/1000-profile inputs follow the shared package and paper-analysis structure. See [Ocean instructions](paper/applications/README_ocean.md) for data preparation, fitting and reporting.

Gaussian process regression with ordinal measurements and potentially available
continuous measurements associated with those inputs. The package provides
fitting, prediction, imputation, diagnostics, continuation, and reporting.

## Code structure

```
eivGP/R/
  helpers.R     # data preparation, kernels, prediction, storage, workflow helpers
  mcmc.R        # posterior sampling, fitting, continuation, sampler diagnostics
  plotting.R    # summaries, comparison tables, figures, report generation
data/ocean_silicate/
  n200/         # frozen 200-profile data and five-fold inputs
  n500/         # updated 500-profile data and five-fold inputs
  n1000/        # raw selected 1000-profile cohort
  prepare_folds.R # standalone fold preparation entrypoint
paper/
  numerical/    # study1 and study2: one .R, .Rmd, and .html each
  applications/ # ocean and adni: one .R, .Rmd, and .html each
  config/       # Ocean dataset/compute parameter templates
  manifests/    # selected numerical runs and frozen Ocean input identities
tools/
  build_package.R
  render_reports.R
```

`eivGP/R/` is the authoritative implementation. The three implementation modules are installed together as one package, following the upstream structure. Metadata and static-analysis declarations remain separate. Analysis scripts call package APIs and select
settings, inputs, and output locations. Competitor adapters call `kergp`, `LVGP`,
and `EzGP`. Reporting functions also live in the package.

## Install

Use R >= 4.1. Install dependencies once, then install this checkout:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "knitr", "patchwork", "posterior",
                   "TruncatedNormal", "digest", "rmarkdown", "testthat"))
# Needed when fitting the competitor methods:
install.packages(c("kergp", "LVGP", "EzGP"))
```

```sh
Rscript tools/build_package.R
R CMD INSTALL eivGP
```

The build helper refreshes the module identity manifest and checks parsing. It
never generates or replaces the package source. Reports require Pandoc, such as
the version included with RStudio.

## Reproduce and read the paper analyses

See [paper/README.md](paper/README.md) for the four workflows and exact commands.
The numerical HTML reports read the selected saved paper exports. The Ocean HTML examples use short-budget results; full-budget application findings remain pending. Rendering does not run MCMC unless `action="run"` is explicitly selected.

Large data, fits, restart checkpoints, and caches stay in the external results
folder configured by `EIVGP_RESULTS_ROOT`. Existing results are not migrated or
rewritten by installation. New code identities require separate run identities;
old caches must not be renamed to force reuse.

## Develop and test

Edit the three package modules, then refresh the manifest and reinstall. API
exports and help pages are maintained in `NAMESPACE` and `man/`.

```r
library(eivGP)
testthat::test_dir("eivGP/tests/testthat", env=new.env(parent=asNamespace("eivGP")),
                   stop_on_failure=TRUE)
```

Version 0.4.0 consolidates the implementation and analysis/reporting workflows.
It retains the 0.3.1 sampler and scientific defaults. Package tests cover both
engines, prediction, diagnostics, continuation, and lossless compact storage.
