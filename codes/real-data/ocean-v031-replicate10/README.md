# Ocean eivGP 0.3.1 replicate-10 bundle

This folder is the paper alignment for the ocean N/P coarsening experiment
under **eivGP 0.3.1**. It does **not** replace
`codes/real-data/01_ocean_representative_figures.R` or
`02_ocean_replicates.R`. Those older scripts still use the Study-1 embedding
comparators (GP-LearnedEmb / CondMean / Gaussian). The numbers in the current
manuscript table come from the runners here.

The frozen 100/300 cohort is **not** duplicated here. Scripts read

```text
codes/real-data/ocean_data/study1_data.csv
codes/real-data/ocean_data/meta.json
```

Keep this folder next to `ocean_data/`, or set `OCEAN_DATA_DIR`. Nested
`|O|` redraw flags stay in `calib_designs/` because they are specific to this
10-redraw design.

Pinned package commit:

```text
a3aa0f697f97d4bdf027507ac3c9b952644affb9
```

Install from the repository root, then confirm `packageVersion("eivGP")` is
`0.3.1`.

## Design (not k-fold)

- Frozen split: 100 train / 300 test in `../ocean_data/study1_data.csv`.
- Six ordered phosphate classes; test methods receive only `(x, c)`, where
  `x = (temperature, salinity)`.
- Nested calibration sizes `|O| in {10, 25, 50}` with 10 redraws of which
  training rows reveal exact phosphate. Flags are in `calib_designs/`.
- `|O| = 0` is a separate fixed-design single fit, not part of the redraw grid.
- This is **not** 3-fold CV.

Full EIV-GP settings: 4 chains, 20,000 iterations, 5,000 burn-in, 600
predictive draws, squared-exponential kernel. `|O| = 0` uses 30,000 iterations
and 7,500 burn-in.

## Layout

```text
codes/real-data/
├── ocean_data/                         # existing GitHub cohort (reused)
│   ├── study1_data.csv
│   └── meta.json
└── ocean-v031-replicate10/
    ├── 01_run_replicate_one.R          # one redraw × one |O|
    ├── 02_aggregate_replicates.R
    ├── 03_run_fixed_single.R           # |O| = 0/10/25/50
    ├── 04_run_published_baselines.R    # UC-GP, LVGP
    ├── 05_run_ezgp_wrapper.R
    ├── 06_run_lm_class.R
    ├── run_replicate_one.sh
    └── calib_designs/
```

## Smoke check

From this folder:

```sh
# One redraw, |O|=10, short MCMC
Rscript 01_run_replicate_one.R \
  ../ocean_data calib_designs outputs/replicate10 1 10 smoke 0 2

# Linear baseline (reads ../ocean_data by default)
Rscript 06_run_lm_class.R

# UC-GP / LVGP (EzGP author predictor is expected to fail on this split)
Rscript 04_run_published_baselines.R
```

Or:

```sh
./run_replicate_one.sh 1 10 smoke 2
```

## Paper runs

```sh
export EIVGP_PACKAGE_COMMIT=a3aa0f697f97d4bdf027507ac3c9b952644affb9
export EIVGP_CHECKPOINT_EVERY=2000

# 10 nested redraws at |O| = 10, 25, 50
for n in 10 25 50; do
  for r in $(seq 1 10); do
    Rscript 01_run_replicate_one.R \
      ../ocean_data calib_designs outputs/replicate10 "$r" "$n" full 0 4
  done
  Rscript 02_aggregate_replicates.R "$n"
done

# Fixed-design EIV-GP, including the uncalibrated |O|=0 row
for n in 0 10 25 50; do
  Rscript 03_run_fixed_single.R \
    "outputs/fixed_o${n}" full 4 "$n"
done

# Comparators
Rscript 04_run_published_baselines.R
Rscript 05_run_ezgp_wrapper.R
Rscript 06_run_lm_class.R
```

Completed EIV redraws write `complete.flag` and are skipped on rerun. Resume
an interrupted MCMC from `fit_v031.rds`.

Set `TARGET_DIAGNOSTICS=1` (7th argument of the redraw runner) only when a
package `diagnose_eivgp()` audit is needed; it is extra computation.

## Outputs used in the manuscript table

| Row | Script | Typical metrics file |
|---|---|---|
| EIV-GP \|O\|=10/25/50 | `01` + `02` | `outputs/replicate10/aggregate/c*/prediction_mean_sd.csv` |
| EIV-GP \|O\|=0 | `03` | `outputs/fixed_o0/full/predictive_metrics.csv` |
| UC-GP, LVGP | `04` | `outputs/published_baselines/baseline_predictive_metrics.csv` |
| EzGP (wrapper) | `05` | `outputs/ezgp_wrapper/predictive_metrics.csv` |
| LM-Class | `06` | `outputs/lm_class/predictive_metrics.csv` |

The author-package EzGP predictor returns negative variances on four test
records whose `(temperature, salinity, c)` also occur in the training set.
Script `05` keeps that author fit and scores the covariance-consistent wrapper.
