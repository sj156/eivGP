# Ocean silicate data

| Folder | Cohort | Included files |
|---|---:|---|
| `n200` | 200 profiles | Cohort, five train/test pairs and test-u masks |
| `n500` | 500 profiles | Cohort, five train/test pairs and test-u masks |
| `n1000` | 1000 profiles | Raw selected cohort and source metadata |

Profiles are cruise–station units. The application uses salinity/log pressure x,
oxygen u, silicate y and an imposed ordinal recording rule c with boundaries
0, 50, 100, 150 and 180 micromol/kg, within pressure 300–1500 dbar.

The 500 cohort is `seed_465000_joint_quintile_pressure`; the 1000 cohort
is `seed_144069_oxygen_interleaved`. These are development cohorts from the earlier
screening exercise; source identifiers and sampling metadata are retained.

## Calibration and folds

The supplied 200/500 five-fold inputs are ready to run. Within each ordinal
class, training calibration selects 20%/50% of profiles with O20 nested in O50.
Selection balances empirical-quantile coverage of salinity and log pressure,
using recorded C/X without exact oxygen or silicate. This is approximate Latin
coverage; finite candidates, tied values and nesting can prevent strict Latin
margins. `calibration_X_coverage.csv` records the achieved coverage, and
`freeze_manifest.json` records seeds and checksums.

To generate a separate input directory using the supplied partitions:

```sh
Rscript tools/ocean-prepare-folds.R ocean500 /external/ocean500-inputs
```

`ocean200` uses the same interface. Preparation retains cohort/fold assignments
and test masks and computes `calib_o20`/`calib_o50`. Select the directory with
`--data-root` for every method/report command. Existing destinations are not
overwritten.

For the raw 1000-profile cohort, prepare folds before fitting:

```sh
Rscript tools/ocean-prepare-folds.R ocean1000
```

The output is `reproduction/ocean/prepared/ocean1000`, automatically selected by
`--dataset=ocean1000`. The 500 fold method forms oxygen quintiles and
pressure-ordered quintets within C; the 1000 method uses oxygen-ordered quintets.
Every quintet contributes one profile per fold. Assignment balances oxygen,
log pressure and salinity moments. Training-calibration selection uses recorded
C/X only. Test masks use seed `202609260 + fold` and nested class-balanced
20%/50% quotas; the 500 O50 design alternates 12/13 profiles per class.

See [Ocean instructions](../../paper/applications/README_ocean.md) for installation
and fitting/report commands.

## Data attribution

Source: [GLODAPv2.2023 merged and adjusted product](https://doi.org/10.25921/zyrq-ht66).
Lauvset et al. (2024), *Earth System Science Data* 16, 2047–2072,
[product paper](https://doi.org/10.5194/essd-16-2047-2024).
