# Study II: one complete paper workload per machine

The formerly separate follow-up experiments are part of the paper design.
There is **no core/follow-up switch**. The two scripts separate methods so they
can run on different machines, not old and new experiments.

| Script | Complete paper workload |
|---|---|
| `study2-eivGP.R` | Original EIV-GP grid; reuse its calibrated four-proxy fits for independent-test FULL scoring; CUT feedback comparison |
| `study2-competitors.R` | Original UC-GP, LVGP, EzGP; all three calibration-aware complete-case versions |
| `study2.R` / `study2.Rmd` | Read saved results; no fitting |

## Fresh replication

Install this checkout on each machine with `Rscript tools/setup_numerical.R`.
Use matching source and locked dependencies. Plan first:

```sh
Rscript paper/numerical/study2-eivGP.R --action=plan --cores=4 --output="$HOME/eivGP-results/study2-eivGP" --data-root="$HOME/eivGP-results/data/study2"
Rscript paper/numerical/study2-competitors.R --action=plan --cores=4 --output="$HOME/eivGP-results/study2-competitors" --data-root="$HOME/eivGP-results/data/study2"
```

Replace `plan` with `run` on the respective machines for complete reproduction.
Both generate the same fixed-seed datasets and calibration masks locally. The
competitor machine needs no EIV-GP fit archives. The EIV-GP machine automatically
creates a FULL map from its newly saved core run; FULL is **scored again**, not
fitted a second time. Missing or ambiguous archives stop scoring rather than
triggering replacement fits. CUT retains the paper's separate response-free
measurement and conditional-response fits.

The original grid and calibrated comparison budgets remain distinct because
those are the paper's evaluation protocols. The calibrated comparisons use four
chains, measurement 100000/20000 total/warmup, 32 outer states, conditional GP
4000/1000 with four selected states per chain, 128 latent draws, 16 GP-response
and 64 competitor-response draws. Original fits use their fixed paper profile.
CUT is secondary evidence but remains included because its results appear in the
appendix. Pilot and integration-sensitivity variants are not default paper jobs.

## Resume, select methods, or use old fits

Repeat the same command to resume compatible saved work. Preserve output roots,
seeds, budgets and source identity. Historical output folders are not relabeled
or automatically migrated to newer code. Do not update or restart a running job
in the same folder; an existing process continues with its already loaded code.

Method selection is optional. `--methods=CUT` runs only the feedback comparison;
`--methods=UC-GP-CC,LVGP-CC,EzGP-CC` runs only the calibrated competitors. These
selections accept `--reps` and `--calib` for incomplete subsets. The original
paper grid is fixed, so subset flags are rejected when original fitting is also
selected. The audited remaining subsets are calibration 50 / reps 15:50 and
calibration 20 / reps 29:50, totaling 174 competitor and 58 CUT results.

For evaluation of existing FULL fits, explicitly supply `--full-fit-root` or
`--full-fit-map` with `--methods=FULL` (or FULL,CUT). This skips original EIV-GP
fitting. A CSV map has `rep,calib,path` and matching fits are validated before
scoring. This is an optional reuse path, not a requirement for fresh replication.
The already recovered 100 historical FULL scores need no rerun for the existing
paper results.

## Saved reporting and output identity

Use the same script and scientific settings with `--action=report`. The notebook
accepts `lane: eivGP`, `competitors`, or `historical`; no phase parameter is needed.
Set its `results_root` to that lane's output and supply the original options in
`options`. Neither reporting entry point can enable fitting.

Keep the existing internal storage layout for resume compatibility: original
EIV-GP results under `core/<run-id>/`, original competitors under
`original-competitors/`, and calibrated comparisons under
`publication/full/reps-.../reference/mechanism/`. These are storage details,
not separate user workflows. Reports contain raw per-dataset scores, summary and
paired CSVs, failure counts and diagnostics. Native competitor exports also
include LaTeX summaries. The historical notebook reports the audited merged core
export once and labels it separately from the calibrated results.

When collecting machines' results, pair raw scores by study, regime, replication,
calibration and method; validate data identities and settings, and reject duplicate
keys. Never average machine-level means or replace one host's protocol with another.
Cross-lane paired comparisons require collecting both lanes' raw scores first.
Retain failures and diagnostic warnings; replication MCSE excludes within-fit error.
No new fitting was used to validate this workflow.

## Study I

`Rscript paper/numerical/study1.R --action=run --profile=paper` already runs its
original study and all three calibrated competitors together. Its complete-case
paper defaults now match the reported study: 12000/4000 measurement iterations,
16 outer states, 128 latent and 64 response draws. No temporary follow-up launcher
is required for fresh reproduction. Existing processes and saved runs are untouched.

## Convergence and failure reporting

New calibrated runs retain the native optimizer status separately from the
measurement diagnostic pass and record the stage of errors. `optimizer` denotes
the native fitting adapter (including its validation/prediction); it is not proof
that every exception was caused by numerical optimization. Nonconvergence is
reported separately from adapter exceptions. Old errors without stage evidence
remain unknown. No numerical budgets or fitting algorithms change.

Normal calibrated summaries additionally write a sibling `*-convergence-report/`
folder. For existing or running experiments, use the independent saved-artifact
reporter instead of restarting the study or invoking its resume/provenance gates:

```sh
EIVGP_R_LIBRARY=/path/to/separate/reporting-library Rscript paper/numerical/convergence-report.R sources.csv /path/to/separate/report-output
```

Install the new package into that separate library; do not replace a library being
used by running jobs. `sources.csv` has `design,root` columns; each root is the
exact directory containing `job_registry.rds`. Use distinct design labels for
Study I regimes and Study II. To collect machine lanes under one design, first
verify identical data, calibration masks, test targets and scoring settings.
Duplicate method/replication keys are rejected. Original-grid exports use a
different protocol and are not pooled with calibrated-comparison scores.

`failure_rates.csv` gives planned, pending, resolved, complete and screened-pass
counts, and separate optimizer, measurement and response-stage failure rates.
Every stage rate includes its assessed denominator and unknown count. Shared
measurement results count as method-job exposures, not independent fits.
`matched_converged_summary.csv` gives paired means, differences, replication MCSE
and matched sample sizes, using only runs passing every applicable check. Unknown
checks are excluded; finite-score matching is metric-specific. Unfiltered scores
remain available. Interpret matched accuracy alongside failure rates because
conditioning on success can select easier datasets. Diagnostics do not prove
convergence. Render `convergence-report.Rmd` with `manifest` and `output` parameters
for the same analysis in HTML. Partial-run reports are snapshots, not final results.

## Scientific ESS continuation (joint EIV-GP fits)

The paper-profile core FULL workflow uses `config$mcmc$ess_min = 100`.
It runs the initial budget, then checks scientific-target bulk and tail ESS
on all retained draws at the fixed five-point test panel. Compatible saved fits
are checked and extended before prediction/recovery scoring. Additional batches
are estimated from the weakest ESS, with a 20% margin and 5,000–20,000 transitions
per chain per check. Sampling stops only when both ESS types strictly exceed the
threshold; R-hat remains reported, not a stopping condition. Persistent separated
chains still require investigation. Set `config$mcmc$ess_min <- NULL` in a custom
configuration to retain fixed-budget behavior.

For ordinary package fits, automatic continuation is opt-in:

```r
fit <- fit_eivgp(X, y, C, engine = "multivariate", latent_dim = 2,
  ess_min = 100,
  ess_panel = list(X = X_eval, C = C_eval, U = U_eval),
  n_iter = 20000, burn = 5000, n_chains = 4)
```

`ess_min = NULL` is the base API default. Without `ess_panel`, the monitored
observable targets use five evenly spaced training rows. Supply integer-coded C
and a meaningful fixed panel for scientific applications. U is optional and adds
physical-surface targets when calibrated. No new warmup or discarded history is
introduced during continuation. `fit$scientific_ess` stores the checks and traces.
For existing fits, `continue_eivgp_until_ess(..., ess_min=100, adaptive=TRUE)`
uses the same logic, including compatible internal Study II saved fits.

With `--mixing-plots=true` (the default), core reporting produces
`scientific_traces_means.{png,pdf}` and
`scientific_traces_uncertainty.{png,pdf}`, each with up to four weak scientific
targets. Display subsampling spans the complete retained history; ESS uses every
draw. This does not change the separate CUT measurement/conditional-GP sampling
budgets or imply that their nested-integration precision has been established.

Do not reinstall the package during an active run. Source updates take effect
for future sessions after `Rscript tools/setup_numerical.R` has been run from the
repository root. A process already running keeps its loaded implementation.


### CUT: descriptive scientific traces, no ESS stopping requirement

CUT remains a fixed-budget comparison method. Neither outer nor conditional
ESS triggers continuation or excludes completed scores. Diagnostic flags remain
visible and do not certify final CUT precision. The report now includes two
separate trace groups: full retained conditional-GP histories for the first
outer state of the first selected replication at each calibration, and the
retained outer integration sequence of conditional posterior means and
predictive second moments. Conditional plots hold the outer U/A/threshold state
fixed; incidental measurement-parameter draws from the GP fit are not substituted
for it. Outer summaries average the available inner GP draws and retain their
finite Monte Carlo error; they are not every measurement-chain transition.
The plots and diagnostics are generated from saved draws, without extra MCMC.
