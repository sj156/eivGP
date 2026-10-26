# Ocean silicate data

| Folder | Cohort | Included files |
|---|---:|---|
| `n200` | 200 profiles | Cohort, five train/test pairs and test-u masks |
| `n500` | 500 profiles | Updated cohort, five train/test pairs and test-u masks |
| `n1000` | 1000 profiles | Raw selected cohort and source metadata |

Profiles are cruise–station units. The application uses salinity/log pressure x,
oxygen u, silicate y and an imposed ordinal recording rule c with boundaries
0, 50, 100, 150 and 180 micromol/kg, within pressure 300–1500 dbar.

The updated 500 cohort is `seed_465000_joint_quintile_pressure`; the 1000 cohort
is `seed_144069_oxygen_interleaved`. These are development cohorts from the earlier
screening exercise; source identifiers and sampling metadata are retained.

## Fold preparation

From the repository root:

```sh
Rscript tools/ocean-prepare-folds.R ocean1000
```

`tools/ocean-prepare-folds.R` performs the input-balanced assignment, nested training calibration
and test-u mask generation entirely in R. It uses `jsonlite` for metadata.
Output defaults to `reproduction/ocean/prepared/ocean1000` and is selected by
`--dataset=ocean1000` in `ocean.R` and `dataset="ocean1000"` in the Rmd. Pass a second argument to choose another
output directory. `ocean500` can also be regenerated with the same command.

R uses its own random-number generator with fixed seeds. New partitions are
reproducible within this R implementation; the supplied 200/500 folds are preserved.

The 500 method forms oxygen quintiles and pressure-ordered quintets within each c;
the 1000 method forms successive oxygen-ordered quintets. Every quintet contributes
one profile per fold. Assignment and within-quintet swaps balance U, log pressure
and salinity moments. The supplied 200/500 calibration files retain their historical
oxygen-boundary anchors. New preparation uses C-stratified, nested empirical-rank
coverage in salinity and log pressure only; no exact U, Y or U-derived cell is
used to choose calibration. Finite candidates, ties and nesting can prevent strict
Latin margins; `calibration_X_coverage.csv` discloses the achieved coverage. Test masks use R seed `202609260 + fold` and class-balanced
20%/50% quotas; the 500 design alternates the 12/13 O50 class quotas.

See [Ocean instructions](../../paper/applications/README_ocean.md) for installation
and fitting/report commands.

## Revised 500-profile calibration

```sh
Rscript tools/ocean-prepare-folds.R ocean500 --calibration-only
```

This preserves supplied cohort/fold assignments and test masks, replacing only
`calib_o20`/`calib_o50` in copied training files. Output defaults to
`reproduction/ocean/prepared/ocean500_xcalibration_v1`; existing destinations are
not overwritten. Use this path with `--data-root` and a new `--output` for all
method/report commands. The source CSVs are unchanged. The application guide
contains the complete command sequence.

## Data attribution

Source: [GLODAPv2.2023 merged and adjusted product](https://doi.org/10.25921/zyrq-ht66).
Lauvset et al. (2024), *Earth System Science Data* 16, 2047–2072,
[product paper](https://doi.org/10.5194/essd-16-2047-2024).
