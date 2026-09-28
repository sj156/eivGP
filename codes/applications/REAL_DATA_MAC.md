# Mac mini

From the repository root, with private cleaned inputs already in `real-data/`:

```sh
Rscript codes/applications/setup_real_data.R  # First setup only
Rscript codes/applications/check_real_data.R
Rscript codes/applications/run_real_data.R adni --cores 12 --plan
caffeinate -i Rscript codes/applications/run_real_data.R adni --cores 12
```

`caffeinate` is optional and prevents idle sleep during the R process. No shell
script is required. Replace 12 with your total core budget. The notebook runs
one fixed five-fold CV; results go to
`real_data_application_outputs/adni/five_fold/`. Rerun the same command to
resume; never launch two coordinators against the same output directory.

To inspect a short pipeline run first:

```sh
EIVGP_SMOKE_TEST=1 Rscript codes/applications/run_real_data.R adni --cores 12
```

Smoke files go to `real_data_application_outputs/adni/smoke/`.
See `REAL_DATA_APPLICATIONS.md` for ocean, report rendering, and method details.
