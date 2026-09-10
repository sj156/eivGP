# ADNI Toledo EIV-GP run bundle

This folder is self-contained for the final EIV-GP runs on the frozen
Toledo-inspired ADNI cohort. It intentionally excludes the earlier ordinary
screening gates.

## Default run

`ADNI_Toledo_EIVGP.Rmd` defaults to **repeat 2, folds 1--3**, four chains,
20,000 transitions per chain (4,000 warm-up), and one 10,000-transition
extension when the lightweight convergence gate fails. The fit uses eivGP
0.3.1, `ident="none"`, the SE/ARD kernel, `u_block_size=8`, and
`cross_every=10`.

Open the Rmd in RStudio and click **Knit**, or run from this directory:

```r
rmarkdown::render("ADNI_Toledo_EIVGP.Rmd")
```

On macOS/Linux the four chains use four forked workers. The package's current
backend falls back to serial execution on Windows.

## Change to repeat 3

At the top of the Rmd, change:

```r
DEFAULT_REPEAT_ID <- 2L
```

to `3L`. Alternatively, do not edit the file and launch with
`ADNI_REPEAT_ID=3` in the environment.

## Quick installation and smoke test

The Rmd installs the tested GitHub revision of `sj156/eivGP` automatically when
the package is absent. `rmarkdown`, `knitr`, and `remotes` should be installed
before rendering:

```r
install.packages(c("rmarkdown", "knitr", "remotes"))
```

To check the complete data -> fit -> checkpoint -> diagnostic plots ->
prediction path without starting production:

```sh
EIVGP_SMOKE_TEST=1 Rscript -e 'rmarkdown::render("ADNI_Toledo_EIVGP.Rmd")'
```

Smoke output is isolated under `outputs/smoke_repeat2/` and cannot be mistaken
for production output.

## Checkpoints and progress

- First checkpoint: iteration 6,000 (4,000 warm-up + 2,000 retained).
- Later checkpoints: every 2,000 transitions through 20,000, and through
  30,000 only when the 20k lightweight gate fails.
- `fit_latest.rds` is atomically overwritten. It contains every retained draw
  and the continuation state, so keeping older copies is unnecessary and would
  consume substantial disk space.
- `PROGRESS.md` and `CURRENT_STATUS.txt` live in each fold output directory.
- Re-running the Rmd resumes a compatible `fit_latest.rds` rather than starting
  that fold over.

The package cannot expose a recoverable checkpoint before a `fit_eivgp()` or
`continue_eivgp()` call returns; therefore iteration 2,000 and 4,000 cannot be
saved while retaining a 4,000-transition warm-up.

## Convergence outputs

The workflow deliberately does not run the expensive full
`diagnose_eivgp()`. It uses the raw diagnostics already computed in the fit and
writes:

- parameter and summary diagnostic CSVs at 20k and, if needed, 30k;
- key-parameter trace plots;
- trace plots for the six missing-U coordinates with the worst raw R-hat;
- an R-hat overview plot with the 1.01 threshold.

The lightweight gate requires finite diagnostics, R-hat <= 1.01, bulk and tail
ESS >= 400, and MCSE/SD <= 0.10 for continuous hyperparameters, loadings,
thresholds, and every missing U. Failure at 30k is reported, not extended again.

## Frozen files

- `data/toledo_adni_cohort_n495.csv`: full n=495 cohort.
- `data/toledo_adni_balanced_repeated_3fold.csv`: frozen repeat 1--3 folds.
- `DATA_MANIFEST.txt`: checksums and design constants.

Do not regenerate or edit the CSVs. The Rmd checks their MD5 hashes before any
fit begins.

