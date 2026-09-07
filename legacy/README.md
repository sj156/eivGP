# Legacy literate chapters

These files preserve the pre-refactor Study I package and analysis chapters.
They are intentionally excluded from `_bookdown.yml` because `litr` extracts
roxygen-led chunks even when ordinary knitr evaluation is disabled. Including
the chapters in the active build would regenerate stale and duplicated package
functions.

The root `.Rmd` book chapters, `_bookdown.yml`, and their payload under
`script/newly-written-and-modified-scripts/` are also historical 0.1.x source.
They must not regenerate the current package. The sole active 0.3.0 build is
`Rscript --vanilla litr/render-package.R`, which uses canonical `codes/`
modules and writes `eivGP/`. The retired 0.2.1 build is recorded in
`litr/legacy/`. Older R code also remains in `script/legacy-code-scripts/`.
