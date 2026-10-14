# 1000-profile raw subset

Candidate: `seed_144069_oxygen_interleaved`, sampling seed 144069.
`cohort_manifest.csv` contains the selected measurements for 1000 distinct
profiles, 250 per ordinal class. Fold assignments, train/test files and masks
are omitted. Source and selection metadata are in `cohort_metadata.json`.

Prepare five folds before using `dataset="ocean1000"` in the Rmd:

```sh
Rscript data/ocean_silicate/prepare_folds.R ocean1000
```

Run from the repository root. Prepared files are written to
`reproduction/ocean/prepared/ocean1000`. Each fold has 800 training and 200 test
profiles, training calibration sizes 0/160/400 and measured test-u sizes 40/100.
See the [data guide](../README.md) for the method and the
[Ocean instructions](../../../paper/applications/README_ocean.md) for fitting.
