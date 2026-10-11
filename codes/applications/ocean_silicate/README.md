# Ocean silicate application

This application models silicate using salinity, log pressure and oxygen. Oxygen is supplied as one of four recorded bands, with exact measurements available for a subset of the training profiles. A parameterized R Markdown document runs selected fits and assembles the analysis into a self-contained HTML report.

## Repository layout

Paths below are relative to the repository root.

| Path | Purpose |
|---|---|
| `replication/02_ocean.Rmd` | Parameter selection, execution and HTML presentation |
| `codes/applications/ocean_silicate/01_setup.R` | Fixed settings and seeds, frozen-input checks, task selection, continuation decisions, resume and merge checks |
| `codes/applications/ocean_silicate/02_fit_models.R` | EIV-GP fitting, continuation, prediction and imputation; native baseline fitting and prediction |
| `codes/applications/ocean_silicate/03_reporting.R` | Execution logging, checked result merging, prediction summaries and paper tables |
| `codes/applications/ocean_silicate/04_diagnostics.R` | U-imputation diagnostics, category means, continuous response functions and measured-profile overlays |
| `data/ocean_silicate/` | Cohort manifest, five training CSVs, five test CSVs and five test-U mask CSVs |
| `docs/replication/02_ocean.html` | Rendered report |
| `reproduction/ocean/` | Locally generated fits, logs, tables and figures |

The Rmd sources the four numbered scripts. Sourcing defines functions without fitting. The marked package calls in `02_fit_models.R` are also displayed in the HTML: `eivGP::fit_eivgp()`, `eivGP::continue_eivgp()`, `eivGP::predict_eivgp()` and `eivGP::impute_eivgp()`. Native baseline calls use the pinned package's internal adapter, `eivGP:::run_published_mixedgp_competitors()`, which invokes kergp, LVGP and EzGP fit/predict functions.

## Environment

Use R with the following modeling packages: eivGP 0.3.1, kergp 0.5.8, LVGP 2.1.6 and EzGP 0.1.0. Also install digest, ggplot2, patchwork, knitr and rmarkdown. HTML rendering requires Pandoc; RStudio normally supplies it. The workflow does not install packages automatically.

Use this fixed eivGP source to obtain the required public functions and baseline adapters:

```r
remotes::install_github(
  "sj156/eivGP",
  ref = "671ce98339807cd3a13347e8468f36acf097d50a",
  subdir = "eivGP", upgrade = "never"
)
rmarkdown::pandoc_available()
```

Prepare the listed dependency versions before rendering. Package version numbers alone do not identify every function implementation. Runtime result identities therefore record package versions, function fingerprints, input/code hashes and numerical settings automatically. Separate environment or checksum CSVs are not required to run the Rmd.

## Run the analysis

Run the examples from an R session at the repository root. Keep the Rmd in `replication/`: it derives the repository root from that location, and resolves input/result paths relative to that root.

The default action is `plan`, which validates the frozen design and displays selected tasks without fitting:

```r
rmarkdown::render("replication/02_ocean.Rmd",
  params = list(profile = "paper", action = "plan"))
```

Run the complete experiment and write the report to its repository location:

```r
dir.create("docs/replication", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(
  "replication/02_ocean.Rmd",
  output_file = "02_ocean.html", output_dir = "docs/replication",
  params = list(profile = "paper", action = "run",
    output_dir = "reproduction/ocean/full", run_id = "full")
)
```

The render argument `output_dir` selects the HTML destination. `params$output_dir` selects the result root containing fits and diagnostic outputs. Clicking Knit with the default parameters produces a plan; set `action` to `run` to fit models.

For a small workflow check, use `profile="quick"`, `action="run"`, `cores=1`, and `modules=c("eiv","prediction")` with a separate result root. Quick fits use fold 1/O0/O50, one chain and 80 iterations; they are not paper results and do not produce the full diagnostics.

## Select and distribute work

| Parameter | Meaning |
|---|---|
| `action` | `plan`: inspect tasks; `run`: execute/resume; `report`: read saved results |
| `profile` | `paper`: full numerical settings; `quick`: small workflow check |
| `folds`, `cells` | Folds 1–5 and calibration cases O0/O20/O50 |
| `modules` | `eiv`, `original_baselines`, `continuation`, `augmented_baselines`, `prediction`, `diagnostics`, `surface_comparison` |
| `augmented_methods`, `placeholders` | Select UC-GP/LVGP/EzGP and expanded-input encoding; primary placeholder is 0 |
| `task_ids` | Optional exact task selection within the other selectors |
| `cores`, `diagnostics_workers` | Fit and diagnostic parallelism |
| `output_dir`, `run_id` | Result root and log identifier; neither changes seeds |
| `import_dirs` | Compatible completed result roots to merge |
| `representative_fold`, `appendix` | Main-figure fold and all-fold supplements |

For example, assign folds 1–3 to one machine:

```r
rmarkdown::render("replication/02_ocean.Rmd", params = list(
  profile = "paper", action = "run", folds = c(1, 2, 3),
  modules = c("eiv", "original_baselines", "continuation",
              "augmented_baselines", "prediction"),
  output_dir = "reproduction/ocean/part_a", run_id = "part_a"))
```

On another machine, use `folds=c(4,5)`, `output_dir="reproduction/ocean/part_b"` and `run_id="part_b"`. Both machines must use the same scripts, frozen data and package implementations. Transfer both complete result roots to one machine, preserving their internal structure and receipts, then merge and report:

```r
dir.create("docs/replication", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("replication/02_ocean.Rmd",
  output_file = "02_ocean.html", output_dir = "docs/replication",
  params = list(profile = "paper", action = "report",
    output_dir = "reproduction/ocean/combined",
    import_dirs = c("reproduction/ocean/part_a", "reproduction/ocean/part_b"),
    modules = c("prediction", "diagnostics", "surface_comparison"),
    run_id = "combined"))
```

Use a fresh combined root. Import validates completed task receipts and identities before copying; conflicting duplicate tasks are rejected. A partial run reports its completed fold count. Full U/m/f displays require all five original baseline groups, all 15 final EIV cells and their diagnostic-selected continuations. Change selectors to add tasks; validated completed tasks are reused. Use a new result root when changing scientific code, inputs or package implementations.

## Design and numerical settings

The frozen cohort contains 200 distinct North Pacific profiles, split into five folds of 160 training and 40 held-out profiles. Salinity and log pressure are fully observed. Oxygen is below 180 micromol/kg and recorded through boundaries 50, 100 and 150. Training cases O0/O20/O50 reveal 0/32/80 exact oxygen values. Nested mixed-test masks reveal 8 or 20 exact values; held-out Y is never used for fitting or tuning.

Each EIV cell initially uses four chains of 20,000 iterations, burn-in 5,000 and `u_block_size=8`. Its current initial diagnostics trigger one additional 20,000-iteration continuation if maximum finite Rhat exceeds 1.01, minimum finite bulk ESS is below 400, or no finite values are available. Passing candidates are skipped. Undefined diagnostics and tail ESS are retained for inspection; remaining failures after the one allowed continuation are reported. Continuation restores saved chains and RNG states.

The primary execution has 50 base tasks plus the selected continuations: 15 initial EIV cells, five original three-method baseline groups, and 30 expanded baseline fits. All 15 cells are continuation candidates. Original baselines use `(x,c)`; expanded baselines use continuous `(x,W)` and categorical `(c,R)`, where R indicates observed oxygen, observed W is standardized using training calibration values, and hidden W takes the primary placeholder 0. Native kernels are used. Class-only and information-matched mixed test predictions are reported separately. Optional placeholders -1 and +1 are available through selectors; primary tables and plots use 0.

## Results and diagnostics

- **Y prediction:** original and expanded baseline comparisons, fold mean/SD, completion counts, raw scores, and compact/full paper tables. Fold SD describes variation across overlapping training sets.
- **U imputation:** bias, RMSE, MAE, 95% coverage and interval width on the same 80 commonly hidden training rows in O20/O50. Figures include O0/O20/O50 and hollow observed-calibration points. O0 uses its unanchored model coordinate and has no physical-scale recovery score.
- **Category means m(x,c):** package category-mean posteriors, original native baseline means and nearby measured profiles; the pressure slice fixes salinity at 34.20.
- **Continuous f(x,u):** package response-function posteriors at salinity 34.25 and pressure 650 dbar, with nearby measured profiles. The combined f/m figure adds baseline category means on its right panel. These noisy measured responses supply data context rather than true functions for surface-recovery RMSE.

Prediction exports are under `prediction_report/`. Diagnostics are under `diagnostics/`, including `tables/paper_u_imputation_main.*`, `tables/paper_u_imputation_class.*`, `prediction/paper_y_prediction_main.*`, `prediction/paper_y_prediction_full.*`, PNG/PDF figures and a lightweight `preview.html`. Paper tables are exported as CSV, Markdown, LaTeX and HTML. Task status, progress and execution logs are under `jobs/<run_id>/`; the HTML includes `sessionInfo()`.

To regenerate the report from completed fits without fitting:

```r
rmarkdown::render("replication/02_ocean.Rmd", params = list(
  profile = "paper", action = "report",
  output_dir = "reproduction/ocean/full", run_id = "report"))
```

For standalone diagnostics, save `list(output_dir="reproduction/ocean/full", report_dir="reproduction/ocean/full/diagnostics")` to a parameter RDS, then run `Rscript codes/applications/ocean_silicate/04_diagnostics.R diagnostics_settings.rds`. This entry does not need Pandoc.

## Data and repository outputs

Data come from the GLODAPv2.2023 merged and adjusted product: DOI [10.25921/zyrq-ht66](https://doi.org/10.25921/zyrq-ht66). Product reference: Lauvset et al., *Earth System Science Data*, 16, 2047–2072, DOI [10.5194/essd-16-2047-2024](https://doi.org/10.5194/essd-16-2047-2024).

The source-region limits are 10–50 degrees N and 140–240 degrees E. Source preparation requires valid salinity and temperature, positive oxygen, and quality flag 2 for oxygen, silicate and nitrate. The cohort further restricts pressure to 300–1500 dbar and oxygen to below 180 micromol/kg. `Y` is silicate, `U` is measured oxygen, and covariates are `G2salinity` and `log_pressure`. Scripts derive the four oxygen classes from U; inherited `c`/`c0` fields are not used. Training calibration masks are `calib_o20`/`calib_o50`; fixed test masks are checked against the prespecified seeds.
