# Paper analyses

Each analysis has one batch `.R` entrypoint, one executable `.Rmd`, and its
rendered `.html`. All reusable computation and reporting comes from `eivGP`.

| Analysis | Batch | Notebook | Read report |
|---|---|---|---|
| Study I | [R](numerical/study1.R) | [Rmd](numerical/study1.Rmd) | [HTML](numerical/study1.html) |
| Study II | [R](numerical/study2.R) | [Rmd](numerical/study2.Rmd) | [HTML](numerical/study2.html) |
| Ocean | [R](applications/ocean.R) | [Rmd](applications/ocean.Rmd) | [HTML](applications/ocean.html) |
| ADNI | [R](applications/adni.R) | [Rmd](applications/adni.Rmd) | [HTML](applications/adni.html) |

## Read the selected paper results

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

## Run new numerical experiments

Each batch script accepts `plan`, `data`, `fit`, or `run`; its default profile is
`paper`. `fit` performs data preparation and our model fitting; `run` also
reports results and runs the competitor packages. Full paper runs are expensive.
Inspect settings with `plan` first:

```sh
Rscript paper/numerical/study1.R plan
Rscript paper/numerical/study1.R run
```

For a quick workflow check, or to set the CPU budget explicitly:

```r
source("paper/numerical/study1.R")
result <- study1_analysis(Sys.getenv("EIVGP_RESULTS_ROOT"), action="run",
                          profile="quick", cores=1L)
```

`quick` skips external competitors and is not an inferential analysis. The
`development` and `paper` profiles use all three competitor packages. Study II
uses the same interface. New runs get code-dependent identities and leave the
selected paper exports intact. See [NUMERICAL_DESIGN.md](../NUMERICAL_DESIGN.md)
for settings. To report a new run, supply its `reporting/combined` directory as
`tables_dir` to the script/notebook; competitor tables can be supplied with
`competitor_dir`. The default report always uses the selected historical export.

## Application workflows

The checked-in application HTML files are plans, not final empirical findings.
Ocean uses the frozen public subset in `data/ocean_silicate`. ADNI requires the
separately supplied collaborator cohort (`toledo_adni_cohort_n495.csv`). Neither
application HTML currently supports a final paper claim; those summaries remain
TODO. Historical application outputs are preserved.

```sh
Rscript paper/applications/ocean.R data/ocean_silicate /external/new-ocean-run plan
Rscript paper/applications/adni.R /external/adni-data /external/new-adni-run plan
```

Replace `plan` with `run` to execute the analysis, or with `report` to process
completed outputs. Render results with explicit notebook parameters:

```r
rmarkdown::render("paper/applications/ocean.Rmd", params=list(
  data_dir=normalizePath("data/ocean_silicate"),
  output_dir="/external/completed-ocean-run", action="report"))
```

The ADNI notebook accepts the same parameters. Use separate output directories
for different package/input identities. Reusable fit storage retains chain
histories and terminal states needed for diagnostics and exact continuation;
existing archives are not automatically converted or deleted.
