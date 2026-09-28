# Real-data application workflows

The application reader entry point is `replication/02_ocean.Rmd`. It provides
exploratory plots, execution, and selected results using eivGP 0.3.1. Its quick
HTML example is in `docs/replication/`. ADNI support is retained for later work;
it is outside the current notebook collection.

## GLODAP silicate/oxygen

`ocean_silicate/` implements the workflow formerly under
`codes/legacy/real-data/ocean_Si_Oxy/scripts/`. Its scientific fitting and
continuation functions are retained. `data/ocean_silicate/` contains the original
200-profile cohort and unchanged five train/test folds.

Quick uses the original script's first-fold O0/O50 demonstration: one chain,
80 iterations, 30 warmup. Paper runs all five folds, O0/O20/O50, and UC-GP/LVGP/EzGP,
with four chains × 20,000 and 5,000 warmup. The recorded selected-chain continuation
adds 20,000 iterations for fold 1/O0, fold 3/O20 and fold 4/O0. Reporting selects
continued results when available and retains original results otherwise.

The older phosphate example is kept in `codes/legacy/ocean_phosphate/`; it is not
the ocean study used by the reader notebook.

## Command-line access

```sh
Rscript codes/applications/setup_real_data.R
Rscript codes/applications/check_real_data.R
OCEAN_QUICK=true Rscript codes/applications/run_real_data.R ocean --cores 4
Rscript codes/applications/run_real_data.R data
```

The notebooks are the primary audience workflow. These optional commands call
the same implementations. Generated outputs remain local; nothing is exported
to Overleaf by rendering a reader notebook.
