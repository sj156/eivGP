# Legacy literate chapters

These files preserve the pre-refactor Study I package and analysis chapters.
They are intentionally excluded from `_bookdown.yml` because `litr` extracts
roxygen-led chunks even when ordinary knitr evaluation is disabled. Including
the chapters in the active build would regenerate stale and duplicated package
functions.

The active package sources are the explicitly ordered root `.Rmd` modules and
the audited scripts under `script/R-scripts/`.
