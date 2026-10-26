# Ocean silicate: 200-profile frozen inputs

This semi-real application predicts silicate y from salinity/log-pressure x and continuous oxygen u or its ordinal record c. Concentrations use micromol/kg. The physical analysis domain is pressure 300–1500 dbar and 0 < u < 180; constructed class boundaries are 0, 50, 100, 150 and 180. Profiles are cruise-station sampling units, not individuals.

## Frozen files

- `cohort_manifest.csv`: 200 distinct profiles, 50 per c level.
- `fold01_train.csv` through `fold05_train.csv`: five training sets, 160 profiles each, with nested calibration masks `calib_o20` and `calib_o50`.
- `fold01_test.csv` through `fold05_test.csv`: five held-out sets, 40 profiles each. Each cohort profile is tested once.
- `fold01_test_U_masks.csv` through `fold05_test_U_masks.csv`: fixed nested test-u masks; 8/20 observed-u test profiles, generated using seed 202609260 + fold. Masks are identical across methods.

Training calibration sizes are 0/32/80. Within each c, the selector uses salinity/log-pressure empirical quantiles for approximate Latin coverage, without exact oxygen or silicate; O20 is nested in O50. Seeds and coverage checks are recorded in `freeze_manifest.json` and `calibration_X_coverage.csv`. The imputation comparison uses the same 80 common-hidden training profiles under O20/O50 (20 per c). Hidden reference U and test Y remain in the files for validation and offline evaluation; the fitting/prediction workflow excludes them as specified in the report.

## Protocols

- **OI:** all training profiles; every test prediction uses (x,c). Native OI comparators ignore calibration u.
- **CC:** native comparators fit only the 32/80 measured-u training profiles and are evaluated on the same 8/20 measured-u test subset. They use continuous (x,u) and categorical c. EIV-GP retains its complete training set and is scored on the identical test subset.
- **MI:** EIV-GP uses all training profiles and calibration u; test prediction routes to (x,u) when u is observed and (x,c) otherwise. All-profile scores include all 40 test profiles.

O0 has an unanchored model coordinate and no physical-scale U recovery score. Marginal distribution plots complement individual interval calibration; marginal spread is not individual uncertainty alone.

## Run

Use `dataset="ocean200"` in [the executable Ocean report](../../../paper/applications/ocean.Rmd). See [Ocean instructions](../../../paper/applications/README_ocean.md). The package validates frozen rows, class counts, calibration masks, nested test masks and train/test separation. Regenerate masks reproducibly in a separate directory with `Rscript tools/ocean-prepare-folds.R ocean200 /external/ocean200-inputs`.

## Data attribution

Source: [GLODAPv2.2023 merged and adjusted product](https://doi.org/10.25921/zyrq-ht66). Cite the dataset and Lauvset et al. (2024), *Earth System Science Data* 16, 2047–2072, [product paper](https://doi.org/10.5194/essd-16-2047-2024). Follow GLODAP fair-data-use and attribution guidance, including cruise-specific references when relevant.
