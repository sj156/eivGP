# GLODAP silicate/oxygen application

Run `replication/02_ocean.Rmd`. The notebook follows the original
`ocean_Si_Oxy` full-five-fold and selected-continuation scripts. This collection
makes their input/output roots explicit so they can be sourced safely from a
notebook; the scientific model calls, seeds, folds, and paper controls are retained.

- `fit.R`: original initial fitting, prediction, imputation, and diagnostics.
- `continue.R`: recorded selected-chain extension and rescoring.
- `workflow.R`: frozen-input checks, source/data cache identity, and orchestration.
- `reporting.R`: combine initial and continued results, including all paper baselines.
- `check_data.py`: independent audit of the bundled cohort and masks.
- `summarize_legacy.py OUTPUT_DIRECTORY`: original detailed Python export of saved results.

The original script quick setting (one chain, 80/30 iterations/warmup; fold 1;
O0 and O50) supplies the compiled vignette. Paper uses all five folds and three
calibration levels, four chains, 20,000/5,000 initial iterations/warmup, the three
published baselines, and the recorded 20,000-step continuation of selected fits.

The unchanged frozen CSVs and source attribution are in `data/ocean_silicate/`.
Fit objects and checkpoints are written beneath the notebook output directory,
not beside these scripts or bundled inputs.
