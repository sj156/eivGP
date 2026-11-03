# Study II: one complete paper workload per machine

The formerly separate follow-up experiments are part of the paper design.
There is **no core/follow-up switch**. The two scripts separate methods so they
can run on different machines, not old and new experiments.

| Script | Complete paper workload |
|---|---|
| `study2-eivGP.R` | Reported fixed-budget EIV-GP grid; independent-test FULL scoring over all 50 datasets; appendix CUT comparison on 28/14 pairs |
| `study2-competitors.R` | Original UC-GP, LVGP, EzGP; all three calibration-aware complete-case versions |
| `study2.R` / `study2.Rmd` | Read saved results; no fitting |

## Fresh replication

Install this checkout on each machine with `Rscript tools/setup_numerical.R`.
Use matching source and locked dependencies. Plan first:

```sh
Rscript paper/numerical/study2-eivGP.R --action=plan --cores=12 --output="$HOME/eivGP-results/study2-eivGP" --data-root="$HOME/eivGP-results/data/study2"
Rscript paper/numerical/study2-competitors.R --action=plan --cores=12 --output="$HOME/eivGP-results/study2-competitors" --data-root="$HOME/eivGP-results/data/study2"
```

Replace `plan` with `run` on the respective machines for complete reproduction.
Both generate the same fixed-seed datasets and calibration masks locally.
The paper scripts default to 12 cores and retain four chains per fit. The
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
selected. Without subset flags, CUT reproduces the reported appendix sample: replications
1:28 at calibration 20 and 1:14 at calibration 50. The exact 42 keys are recorded
in `../manifests/study2/feedback-pairs.csv`. FULL rescoring and the calibrated
competitors still use all 50 datasets per calibration. Explicit subset flags
request a separate subset; extending CUT to 50 changes the appendix experiment.

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
export once and labels it separately from the calibrated results. Historical
CSV-only reporting defaults to `mixing_plots: false`: combined run-level
diagnostic screens cannot supply chain traces. The method-lane reports produce
the saved scientific-target trace figures; explicitly enabling mixing plots in a
historical view requires detailed compatible diagnostic inputs.

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

## Fixed primary budget and optional continuation

The paper lane uses four chains with 20,000 total iterations and 5,000 warmup,
without automatic ESS continuation. Diagnostics and unfavorable fits remain in
the summaries. This restores the budget underlying the reported primary grid;
the older script default that extended fits until ESS >100 produces a different
experiment. Optional continuation remains available through the package API and
must use a separate output root and be labeled as a new analysis.

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


## Compact CUT storage and checkpoint cleanup

New Study II calibrated runs default to `--retain-draws=false` in both lane
scripts. The response-free measurement sampler computes chain-aware diagnostics
and exact imputation means/intervals from every post-warmup draw, then saves only
those summaries, selected outer states, terminal states/RNG, provenance and
verification fingerprints. Imputation scores use these summaries without changing
their definition. The four completed measurement checkpoints are removed only
after the final artifact is saved, read back and verified against all four chains.
Incomplete, locked or mismatched checkpoints remain recoverable. Cached completed
jobs retry cleanup during `run`; `report` does not clean or fit.

Conditional GP histories are temporary working files. Once diagnostics, requested
scientific trace assets and selected integration states have been published, the
verified full conditional fit and redundant restart files are retired. Selected
states and diagnostics remain available for prediction/reports. Use
`--retain-draws=true` from the start to retain full measurement and conditional GP
histories for arbitrary later diagnostics; it cannot reconstruct previously
removed traces. Completed redundant checkpoints are retired with either policy.
The original joint EIV-GP fit archives are still needed for scientific ESS
continuation and the independent FULL evaluation; this CUT retention option does
not remove them or external supplied archives.

For old outputs, use the maintenance script on a **stopped** reference directory
containing both `checkpoints/` and `mechanism/`. It verifies the saved output
rather than restarting sampling, and does not alter historical protocol identity:

```sh
Rscript tools/study2-checkpoints.R --output=/path/to/reference --audit-dir=/local/audits/cut-check-01
Rscript tools/study2-checkpoints.R --output=/path/to/reference --audit-dir=/local/audits/cut-clean-01 --dry-run=false --compact-measurement=true
```

The first command only audits. The second removes verified redundant GP/CUT
measurement checkpoints and, after successful measurement retirement, replaces
the corresponding full measurement trace with compact summaries. Audits must use
fresh local directories outside the results tree. Locks are never overridden.
Historical conditional full fits without the new selected-state verification
metadata remain unchanged. No archive is modified by installing this revision.
New source identity still requires a new run root; maintenance does not relax
cache compatibility or migrate unfinished runs. Install the updated checkout
before using these commands; do not replace a library used by an active process.

## Manuscript selection and reproducibility

The paper profile retains every completed reported setting, including the
category-sufficient boundary and logistic-error pressure test. It omits the older
response-free PI-GP/CC-GP component ablations from new primary jobs; calibrated
literature competitors and the separate Study II CUT comparison are retained.
Removing a duplicated table does not remove its experiment.

The published evaluation budgets, data/calibration seeds and four-chain design
are fixed. The lane scripts reject altered scoring/sampling budgets; use package
APIs and a separate root for sensitivity work. Plan and report do not fit models.
Source changes require new fitting roots; archived results keep their actual
provenance. Install this checkout before a new session.

Reproducing a saved score from the same retained draws should agree to numerical
tolerance. Fresh MCMC fits can additionally differ through Monte Carlo error,
chain exploration and numerical libraries; locked dependencies and seeds do not
justify a floating-point-only guarantee. Compare replication-level paired scores,
diagnostics and numerical precision, not rounded table entries alone. No complete
fresh 50-dataset rerun was used to certify this code revision. Existing HTML files
are historical snapshots; render the matching Rmd against the selected results
for an updated report. FULL-CUT details remain appendix evidence and should have
only a brief interpretation in the main text.
