# 500-profile subset

Candidate: `seed_465000_joint_quintile_pressure`, sampling seed 465000.
The cohort contains 125 profiles per ordinal class. Each fold uses 400 training
and 100 test profiles. Training calibration sizes are 0/80/200; measured test-u
sizes are 20/50. The O50 test mask alternates 12/13 profiles per class.

`cohort_manifest.csv`, the five train/test pairs and five test-mask files retain
the supplied values and assignments. `freeze_manifest.json` records construction.
The same 200 hidden training profiles per fold support the O20/O50 imputation comparison.

Use `dataset="ocean500"` in [ocean.Rmd](../../../paper/applications/ocean.Rmd).
See the [data guide](../README.md) for the fold method and regeneration command.
