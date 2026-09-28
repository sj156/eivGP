# Ocean silicate application

This semi-synthetic application predicts silicate from salinity, log pressure, and dissolved oxygen information. It uses 200 North Pacific GLODAPv2.2023 observations, with one observation per station profile, pressure 300–1500 dbar and oxygen below 180 µmol kg⁻¹. Oxygen classes have boundaries at 50, 100, and 150 µmol kg⁻¹. Five fixed folds each contain 160 training and 40 test observations. EIV-GP is fitted with 0%, 20%, or 50% exact training oxygen and compared with UC-GP, LVGP, and EzGP.

The data and split files are included. Oxygen classes and exact-observation masks are constructed for the experiment; the source oxygen measurements remain in the files.

## Files

- `data/cohort_manifest.csv`: the 200-observation source subset, including original measurements, station identifiers, and the frozen fold assignment.
- `data/fold01_train.csv`–`fold05_train.csv` and matching `*_test.csv`: the fixed train/test sets. Training files include the nested `calib_o20` and `calib_o50` masks.
- `scripts/00_check_data.py`: checks cohort sizes, oxygen classes, train/test separation, and calibration-mask counts.
- `scripts/01_run_full5.R`: fits the three baselines and all 15 EIV-GP models using the installed `eivGP` package; saves predictions, scores, imputation, and MCMC diagnostics.
- `scripts/02_continue_selected_plus20k.R`: uses `continue_eivgp()` to extend fits selected by training-chain diagnostics, then recomputes predictions and scores.
- `scripts/03_summarize_and_plot.py`: combines original and continued results, checks held-out predictions, produces summary tables and diagnostic plots.

Generated results and fit objects are not included in this upload folder.

## Run

Install `eivGP` version 0.3.1 and its package dependencies in R. The scripts call `fit_eivgp()`, `predict_eivgp()`, `impute_eivgp()`, and `continue_eivgp()` directly; the three comparators and scoring helpers are also called from the package. No package source is copied or changed here. The initial fit uses four chains, 20,000 iterations per chain, 5,000 burn-in iterations, and `u_block_size = 8`. The second script adds 20,000 iterations per chain to selected fits.

From this folder, run:

~~~sh
python3 scripts/00_check_data.py
Rscript scripts/01_run_full5.R
Rscript scripts/02_continue_selected_plus20k.R
python3 scripts/03_summarize_and_plot.py
~~~

The Python summary script uses NumPy, pandas, and Matplotlib. Model fits are saved under `full5_class4/`; continuations and summaries are saved under `continuation_plus20k/`. In the completed run, fold 1 at 0%, fold 3 at 20%, and fold 4 at 0% received chain extensions.

Class-only prediction uses salinity, log pressure, and oxygen class. Mixed EIV-GP prediction additionally uses exact oxygen for the designated 20% or 50% of test observations. The comparator implementations use class-only inputs. At 0% calibration, oxygen imputation is reported on the package model scale, not as a raw-concentration error.

## Data source

The subset comes from the [GLODAPv2.2023 merged and adjusted data product](https://doi.org/10.25921/zyrq-ht66). The [GLODAPv2.2023 product paper](https://doi.org/10.5194/essd-16-2047-2024) states that the data are available without restrictions and asks users to follow its fair-data-use and attribution principles. Cite the product paper and dataset when using this subset. For work centered on individual cruises, consult the GLODAP cruise summary table for the contributing cruise DOIs and publications.

Lauvset, S. K., et al. (2024). “The annual update GLODAPv2.2023: the global interior ocean biogeochemical data product.” *Earth System Science Data*, 16, 2047–2072. https://doi.org/10.5194/essd-16-2047-2024.
