# Ocean silicate application

This semi-synthetic application predicts silicate from salinity, log pressure, and dissolved oxygen information. It uses 200 North Pacific GLODAPv2.2023 observations, with one observation per station profile, pressure 300–1500 dbar and oxygen below 180 µmol kg⁻¹. Oxygen classes have boundaries at 50, 100, and 150 µmol kg⁻¹. Five fixed folds each contain 160 training and 40 test observations. EIV-GP is fitted with 0%, 20%, or 50% exact training oxygen and compared with UC-GP, LVGP, and EzGP.

The data and split files are included. Oxygen classes and exact-observation masks are constructed for the experiment; the source oxygen measurements remain in the files.

## Reproduce

Run `replication/02_ocean.Rmd` from the repository root. Its quick and paper
profiles call the scientific code in `codes/applications/ocean_silicate/`.
The CSVs here are unchanged from the original `ocean_Si_Oxy/data/` bundle.
`cohort_manifest.csv` stores the frozen cohort; `fold01_train.csv` through
`fold05_test.csv` store the five fixed splits and nested calibration masks.

The full profile runs four chains, 20,000 initial iterations, 5,000 warmup,
and the recorded selected 20,000-iteration continuations. Generated fits stay
under `reproduction/ocean/`, separate from these inputs.

Class-only prediction uses salinity, log pressure, and oxygen class. Mixed EIV-GP prediction additionally uses exact oxygen for the designated 20% or 50% of test observations. The comparator implementations use class-only inputs. At 0% calibration, oxygen imputation is reported on the package model scale, not as a raw-concentration error.

## Data source

The subset comes from the [GLODAPv2.2023 merged and adjusted data product](https://doi.org/10.25921/zyrq-ht66). The [GLODAPv2.2023 product paper](https://doi.org/10.5194/essd-16-2047-2024) states that the data are available without restrictions and asks users to follow its fair-data-use and attribution principles. Cite the product paper and dataset when using this subset. For work centered on individual cruises, consult the GLODAP cruise summary table for the contributing cruise DOIs and publications.

Lauvset, S. K., et al. (2024). “The annual update GLODAPv2.2023: the global interior ocean biogeochemical data product.” *Earth System Science Data*, 16, 2047–2072. https://doi.org/10.5194/essd-16-2047-2024.
