# Study I: complete paper workloads on separate machines

The calibrated comparisons are part of the main replication workflow. No separate
follow-up launch is needed. The combined `study1.R` remains compatible, but use
these scripts to split work across machines:

| Script | Workload |
|---|---|
| `study1-eivGP.R` | Original Study I eivGP paper grid, with external competitors disabled |
| `study1-competitors.R` | Original UC-GP, LVGP and EzGP plus their calibrated complete-case versions |
| `study1.Rmd` | Saved-result reporting only; fitting actions are rejected |

Install matching source and locked dependencies before starting **new** runs.
Do not replace a library used by running jobs or change their output folders.

```sh
Rscript tools/setup_numerical.R
Rscript paper/numerical/study1-eivGP.R --action=plan --cores=4 --output="$HOME/eivGP-results/study1-eivGP" --data-root="$HOME/eivGP-results/data/study1"
Rscript paper/numerical/study1-competitors.R --action=plan --cores=4 --output="$HOME/eivGP-results/study1-competitors" --data-root="$HOME/eivGP-results/data/study1"
```

Replace `plan` with `run` on the respective machines. The scripts default to plan;
no fitting occurs without an explicit run action. The competitor machine does
not require eivGP response fits. Both machines generate the fixed-seed data from
the same protocol. No new fitting was used to validate these scripts.

Original competitors use all `(X,C,Y)` records and ignore calibrated `U`.
Complete-case competitors use only calibrated `(X,U,C,Y)` records for the response
model, then integrate response-free measurement-model draws at test time.
This is a calibrated adaptation, with its own measurement diagnostics.

The complete-case paper design has 50 replications, calibration sizes 20 and 50,
and both `eta=0` and `eta=1` regimes: 600 method evaluations. Its fixed defaults
are four measurement chains, 12000/4000 total/warmup iterations, 16 outer states,
128 latent draws and 64 response draws. Original fits retain their original paper
profile. Study I does not acquire Study II's CUT ablation or FULL rescoring step.

Optional `--methods=UC-GP-CC,LVGP-CC,EzGP-CC` selects only calibrated competitors;
then `--reps`, `--calib` and `--eta` can select incomplete subsets. Original-grid
fitting uses its fixed paper grid and rejects these subset flags. Repeating a
command resumes compatible saved work; source/protocol mismatches require a new
root, not relabeling old caches. The internal `core/`, `original-competitors/` and
`complete-case/eta*_balanced/` directories preserve the established layout.

## Failure rates and matched diagnostic-screened accuracy

Study I and Study II use the same `report_numerical_convergence()` implementation.
New runs retain optimizer statuses and measurement diagnostic passes separately,
with error-stage metadata. The native fitting-adapter stage includes validation
and preliminary prediction; an exception there is not automatically numerical
nonconvergence. Legacy errors without stage evidence remain unknown.

Normal complete-case summaries also write sibling `*-convergence-report/` folders.
For running or historical jobs, use the independent reader, which avoids the
experiment's fitting/resume gates and writes outside its input roots:

```sh
EIVGP_R_LIBRARY=/path/to/separate/reporting-library Rscript paper/numerical/convergence-report.R sources.csv /path/to/separate/report-output
```

Install the updated package in that **separate** reporting library. The CSV has
`design,root` columns. Each root is an exact directory containing
`job_registry.rds`, e.g. `complete-case/eta0_balanced` or `eta1_balanced`. Assign a
different design label to each regime. Never pool regimes or original-grid and
calibrated evaluation protocols. Only reuse a design label across machine roots
after verifying identical data, masks, targets and scoring settings; duplicate
method/replication keys are rejected.

Outputs include job statuses, stage-specific failure rates with assessed and
unknown counts, pending jobs, unreadable artifacts, unfiltered scores, and paired
accuracy restricted to matched replications passing every applicable screen.
Every paired metric includes its finite matched count, both means, difference,
and replication MCSE. Shared measurement failures are method-job exposures, not
independent failures. Report failure rates beside screened accuracy to expose
selection effects; diagnostic passes do not prove convergence. Empty matched
sets remain empty. No scores are substituted for failed jobs.

Render `convergence-report.Rmd` with `manifest` and `output` parameters for the
same HTML report. `study1.Rmd` accepts `lane: eivGP`, `competitors`, or
`historical`; `options` cannot override its report-only action. Historical views
are labeled snapshots, not the final merged manuscript results.

## EIV-GP report and trace plots

`study1-eivGP.R --action=run` now finishes by preparing the saved-result report
and rendering the existing `study1.Rmd` to `study1.html`. The HTML and CSVs are
written to `core/study1-report/`, or the explicit `--report-output` directory.
The report covers prediction, observed-input means, latent surfaces, imputation,
supplementary controls, completeness, and scientific-target mixing. It retains
all completed scores. Scientific-target bulk and tail ESS must both be **greater
than 100**; R-hat and MCSE remain supporting diagnostics. Historical overall gate
flags are preserved in the source but do not screen the reported accuracy.

To repair reporting after a completed simulation, use `report`, not `run`:

```sh
Rscript paper/numerical/study1-eivGP.R --action=report \
  --run-dir=/path/to/core/study1-publication-v2-RUNID \
  --data-root=/path/to/data/study1 \
  --report-output=/path/to/local/study1-report
```

This route reads the selected run without fitting, resuming, or rewriting it.
An ambiguous parent directory requires an explicit `--run-dir`. It does not
require the historical configuration to match current fitting code. Provenance
records the original saved inputs and the separate reporting implementation.

Three compact trace figures replace diagnostic scatterplots in the EIV-GP lane:
typical mean targets, weakest mean targets, and weakest uncertainty targets.
Graphics run serially, avoiding macOS Objective-C/fork crashes. New fits retain
small chain-resolved scientific-target traces for every replication; full-window
ESS is computed before display thinning. For older archives, selected traces
are reconstructed from retained fits and checksum-verified frozen inputs using
the original target-integration settings. Missing traces are explicitly listed.
The chosen fits and scientific targets are recorded in `trace_selection.csv`.

`study1.Rmd` now defaults to `lane: eivGP`. It accepts `run_dir`, `data_root` and
`report_output`, or an already prepared `report_rds`. Explicit `lane: historical`
and `lane: competitors` remain available. `--mixing-plots=false` skips trace
reconstruction while retaining the numerical report. Rendering requires
`rmarkdown` and Pandoc. No new MCMC is performed during rendering.

### Scientific-target continuation

Production Study I fitting continues saved chains in 5,000-transition increments
until every monitored scientific target has **both bulk and tail ESS >100**.
All original retained draws count. R-hat and traces remain supporting diagnostics;
this precision criterion is not a claim of convergence. Quick runs remain smoke tests.
Final fits retain terminal RNG states. Checkpoint files are removed only after the
complete expected results, final fit histories, and reporting outputs are verified.

To repair an existing run without rerunning acceptable fits or competitors:

```r
eivGP::continue_study1_results(run_dir, data_root = frozen_study1_data,
                              chunk_iter = 5000L, n_cores = 4L)
```

Then use the normal `--action=report --run-dir=...` command to regenerate
`study1.html`. Rendering never continues chains. A single compatible fit can use
`continue_eivgp_until_ess(fit, X, C, U, checkpoint_file = "continuation.rds")`.
The default has no extension-count limit; an explicit `max_extensions` is an
operational limit and returns `success=FALSE` if the requirement is unmet.
