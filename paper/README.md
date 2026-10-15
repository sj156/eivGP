# Paper analyses

Numerical experiments now have a [fresh-checkout workflow](numerical/README.md)
for macOS and Linux: `Rscript tools/setup_numerical.R`, then
`Rscript paper/numerical/study2.R --component=mechanism --action=plan`. It generates synthetic
data and needs no historical exports. The package supplies model fitting;
the numerical runners define the study-specific comparisons and scoring.

The original study/application notebooks and their rendered HTML remain below
as historical report readers.

| Analysis | Batch | Notebook | Read report |
|---|---|---|---|
| Study I | [R](numerical/study1.R) | [Rmd](numerical/study1.Rmd) | [HTML](numerical/study1.html) |
| Study II | [R](numerical/study2.R) | [Rmd](numerical/study2.Rmd) | [HTML](numerical/study2.html) |
| Ocean | [R](applications/ocean.R) | [Rmd](applications/ocean.Rmd) | [HTML](applications/ocean.html) |
| ADNI | [R](applications/adni.R) | [Rmd](applications/adni.Rmd) | [HTML](applications/adni.html) |

## Read the selected historical results (optional)

From the repository root, after installing the package:

```sh
export EIVGP_RESULTS_ROOT="/path/to/ordinal_mixed_input_GP"
Rscript paper/numerical/study1.R report
Rscript paper/numerical/study2.R report
Rscript tools/render_reports.R
```

The renderer reads the numerical exports and renders application plans. It
neither downloads full MCMC objects nor reruns fits. The source result files
remain unchanged. HTML outputs are written alongside their Rmd sources.

The selected sources, relative to the external results folder, are:

- Study I: `study1/runs/recovery-publication-mini/study1` (300 retained original fits).
- Study II: `study2/runs/recovery-publication-linux/study2` (227 original and 73 recovered fits).

[Study I manifest](manifests/study1/manifest.json) and
[Study II manifest](manifests/study2/manifest.json) record source paths, hashes,
and selection rules; their `fits-used.csv` files identify every selected fit.
Read the merged exports once, without appending original fits again. Scores
retain diagnostic warnings. Replication MCSEs and MCMC diagnostics describe
different sources of uncertainty.

## Run numerical experiments from source

Use [numerical/README.md](numerical/README.md) for the dependency lock,
fixed-seed data generation, targeted FULL/CUT/competitor comparisons, resuming
incomplete jobs, and pending manuscript displays. The two original grids use:

```sh
Rscript paper/numerical/study1.R --action=plan --profile=paper --cores=4
Rscript paper/numerical/study2.R --action=plan --profile=paper --cores=4
```

Replace `plan` by `data` to generate synthetic inputs, or by `run` to generate,
fit and aggregate. All default paths are checkout-relative under `results/`.
The legacy positional report API remains available only for explicitly
supplied historical exports; it is not needed for fresh reproduction.

## Application workflows

Ocean uses the same three package implementation modules as the numerical studies and supports 200, 500 and 1000 profiles. See [Ocean instructions](applications/README_ocean.md) for installation, dataset parameters, full-budget execution and saved-report reconstruction. Its `.R` entry calls the same package API as the Rmd.

```r
rmarkdown::render("paper/applications/ocean.Rmd",params=list(
  dataset="ocean500",profile="paper",action="plan"))
```

Set action to run with a fresh external output_dir for fitting, or report with a completed root for reconstruction. Full-budget runs are expensive; pilot and quick budgets only verify the covered workflow. Only ocean.html is versioned as the short-budget example report. Prepare the 1000-profile folds separately using data/ocean_silicate/prepare_folds.R; the Ocean instructions give the command.

ADNI continues to use its existing workflow and separately supplied collaborator cohort:

```sh
Rscript paper/applications/adni.R /external/adni-data /external/new-adni-run plan
```

Use separate output roots for changed package/input identities. Fits retain chain histories and terminal states for diagnostics and continuation. Existing archives are not automatically migrated or deleted.

Targeted numerical additions (complete-case competitors and matched response-free CUT), explicit command-line reruns, checkpoints and scoring are documented in [numerical/README.md](numerical/README.md). The existing study entrypoints also accept the new `--name=value` CLI; their report functions remain backward compatible.
