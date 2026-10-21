# Paper replication

This folder contains the analyses needed for the paper. For Study II, use the
[separate EIV-GP and competitor workflows](numerical/study2-README.md); their
outputs can be produced on different machines and reporting reads saved results.
Study I and application workflows are unchanged.

| Analysis | Batch | Notebook | Read report |
|---|---|---|---|
| Study I | [R](numerical/study1.R) | [Rmd](numerical/study1.Rmd) | [HTML](numerical/study1.html) |
| Study II | [EIV-GP](numerical/study2-eivGP.R), [competitors](numerical/study2-competitors.R), [report](numerical/study2.R) | [Rmd](numerical/study2.Rmd) | [HTML](numerical/study2.html) |
| Ocean | [R](applications/ocean.R) | [Rmd](applications/ocean.Rmd) | [HTML](applications/ocean.html) |
| ADNI | [R](applications/adni.R) | [Rmd](applications/adni.Rmd) | [HTML](applications/adni.html) |

## Read the selected historical results (optional)

From the repository root, after installing the package:

```sh
export EIVGP_RESULTS_ROOT="/path/to/ordinal_mixed_input_GP"
Rscript paper/numerical/study1.R report
Rscript paper/numerical/study2.R
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

See [the numerical guide](numerical/README.md) for Study I and
[Study II reproduction](numerical/study2-README.md) for independent machine
execution, remaining follow-up jobs, output locations and reporting. Default
launches plan the workload; execution requires `--action=run`.

## Application workflows

Ocean uses the same three package implementation modules as the numerical studies.
Run experiments through its `.R` entrypoint, then render `.Rmd` to read saved
results and generate HTML. See [Ocean instructions](applications/README_ocean.md)
for installation, data preparation, execution and report parameters.

```sh
Rscript paper/applications/ocean.R --action=plan --profile=paper --dataset=ocean500 --cores=4
# Use action=fit, then action=analyze; action=continue explicitly extends flagged chains.
```

The default output is `results/ocean/<dataset>/<profile>/`; `--output` overrides
it. The Rmd defaults to `action="report"` and does not fit models. Prepare the
1000-profile folds separately with `data/ocean_silicate/prepare_folds.R`.
Only `ocean.html` is versioned as the short-chain example report.

ADNI continues to use its existing workflow and separately supplied collaborator cohort:

```sh
Rscript paper/applications/adni.R /external/adni-data /external/new-adni-run plan
```

Use separate output roots for changed package/input identities. Fits retain chain histories and terminal states for diagnostics and continuation. Existing archives are not automatically migrated or deleted.

Study II workflow selection and reporting are documented in [its reproduction guide](numerical/study2-README.md).
