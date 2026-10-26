# ADNI application

Use the same three method-oriented entrypoints as Ocean:

All executable scripts are directly in `tools/`. The method entrypoints load
`tools/adni-tools.R`, which contains all ADNI-specific functions, including
fold preparation, fitting orchestration, evaluation, and reporting. The Rmd
report, HTML, and this guide remain under `paper/applications/`.
The installed package supplies general numerical functions; it contains no ADNI
workflow functions. Changes to ADNI scripts do not require reinstalling eivGP.
The former package exports `adni_workflow_api`, `run_adni_application`, and
`run_adni_analysis` are now application functions available after sourcing
`tools/adni.R`. Install the refactored package once to remove the
old package-owned implementation.

```sh
Rscript tools/adni-eivGP.R --plan=true --profile=paper --cores=16
Rscript tools/adni-eivGP.R --profile=paper --cores=16
Rscript tools/adni-competitors.R --profile=paper --cores=16
Rscript tools/adni-report.R --profile=paper --cores=16
```

Run these sequentially against the same input and output directories. The first
method fits compatible missing EIV tasks and immediately scores them. The second
owns original ordinal-input and complete-case UC-GP, LVGP and EzGP tasks and their
scores. The report entry assembles the saved analyses through `adni.Rmd` and
writes a self-contained `adni.html`. Reporting never fits, retries an optimizer
or extends MCMC. `--reuse-only=true` on either method also prohibits fitting.
`--plan=true` lists work without starting it.

The default source directory is `data/adni`; these paper commands save
results in `results/adni`. Supply `--profile=paper` explicitly because the
ADNI CLI defaults to `pilot`, which saves separately in `results/adni_pilot`.
Choose `--cores` for the available machine; it does not alter chain counts.
Fold CSVs and protocols sit directly in the dataset folder, with filenames such
as `fold01_O20_train.csv`; no `inputs/` or per-fold subdirectories are needed.
Both participant inputs and empirical outputs stay out of Git.
A fresh checkout intentionally contains no ADNI records. See the
[data contract](../../data/adni/README.md) before running.

## Where ADNI results are saved

Paper results default to `results/adni/`; `--dataset=3C` uses `results/3C/`.
Other profiles use `results/<dataset>_<profile>/`. `--output` selects another
root. Like Ocean, new empty result roots use task folders directly under it:

```text
results/adni/
  eiv_fold01_O0/             # fit, final predictions, scores and diagnostics
  eiv_fold01_O20/
  eiv_fold01_O50/
  baseline_fold01/           # original UC-GP, LVGP and EzGP
  cc_fold01_O20_LVGP/        # one complete-case competitor task
  data/                     # verified frozen input snapshot
  logs/                     # execution status, worker logs and core events
  reports/eivGP/index.html
  reports/competitors/index.html
  reports/prediction_report/
  reports/diagnostics/
  reports/paper/             # PDF figures, CSV/HTML tables and artifact manifest
  _meta/                    # manifest, source identity and task receipts
  adni.html
```

Other folds/calibrations follow the same names. Earlier evaluation stages stay
inside each EIV task, e.g. `eiv_fold01_O20/stages/S1000/`; final-stage files sit
beside `fit.rds`. Stage labels denote requested budgets, while saved convergence
tables record actual iterations, including optional scientific-ESS continuation.
The marker `_meta/storage-layout.rds` records `adni-flat-v1` and planned stages.

Existing unmarked output roots keep their original paths. They are never moved
or automatically relabeled. The same path interface reads legacy and flat
results. Old default paper folders such as `results/adni_paper/` must be supplied
explicitly with `--output`; no automatic discovery or migration runs. Source
and numerical identity checks still apply to resume operations.

Use the same dataset, profile and input/output flags on all three entrypoints.
The three-item cohort (MEMORY3, MEMORY8 and DIVATT3) uses `data/adni_3C/`:

```sh
Rscript tools/adni-eivGP.R --profile=paper --dataset=3C --cores=16
Rscript tools/adni-competitors.R --profile=paper --dataset=3C --cores=16
Rscript tools/adni-report.R --profile=paper --dataset=3C --cores=16
```

The default nine-item and three-item cohorts have different eligible populations;
they are not a paired comparison of C sets. All participant inputs, snapshots,
fits and empirical reports remain excluded from Git.

## Rebuild the ADNI report without fitting

```sh
Rscript tools/adni-eivGP.R --reuse-only=true --profile=paper --cores=16
Rscript tools/adni-competitors.R --reuse-only=true --profile=paper --cores=16
Rscript tools/adni-report.R --profile=paper --cores=16
```

These commands read compatible saved results. They do not fit missing models,
retry optimizers or extend chains. Workflow completion alone does not establish
convergence.

## Profiles and explicit numerical settings

| Profile | Chains | Total iterations per chain | Burn | Evaluation stages | Native competitor controls |
|---|---:|---:|---:|---|---|
| pilot | 4 | 1500 | 500 | 1000,1500 | Full public paper controls |
| paper | 4 | 50000 | 5000 | 50000 | Full public paper controls |
| quick | 1 | 80 | 30 | 80; fold1 | Reduced workflow-check controls |

ADNI pilot is a five-fold scientific screening budget. It deliberately differs
from Ocean's 80-iteration pilot. The two ADNI stages continue the same chain
states, with one burn. Defaults are HMC at O0/O20 and ESS at O50. Each pilot/paper
score uses 250 selected posterior states per chain. The default task order is
O50, O20, O0, matching the screening workflow. `--cores` caps active workers
through the unmodified public shared CPU pool; resident fit limits remain those
of the public Ocean executor. CPU allocation does not change seeds or chains.

Parameters can be supplied through `--config=FILE` or explicit flags. For example:

```sh
Rscript tools/adni-eivGP.R --profile=paper --cores=16 --iter=20000 --burn=2000 --stages=10000,20000 --output=results/adni_20k
Rscript tools/adni-competitors.R --profile=paper --cores=16 --iter=20000 --burn=2000 --stages=10000,20000 --output=results/adni_20k
Rscript tools/adni-report.R --profile=paper --output=results/adni_20k
```

`--chains`, `--draws-per-chain`, `--first-checkpoint`, `--checkpoint-step`,
`--folds=1:5`, `--cells=O0,O20,O50`, `--sampler-o0=HMC`, `--sampler-o20=HMC` and
`--sampler-o50=ESS` are also available. The first checkpoint must follow burn.
Saved fitting identities freeze numerical and input settings; use a new output
root for a different experiment. Changing cores alone is permitted. The data
protocol determines the number/names of X, U and C coordinates, ordinal levels,
cohort size, fold membership and observation masks; no case-specific N or q is
required in the method scripts.

## Local data and frozen sampling

The current application metadata example describes memory with two imaging
coordinates, nine native ECog items, age and education, within women. A different
compatible item set uses the same methods with its own metadata and prepared
inputs. Changing Y/U definitions requires an explicit scientific design; the
software does not choose a candidate by prediction error.

```sh
Rscript tools/adni-eivGP.R --profile=paper --data-root=data/adni_other --output=results/adni_other_paper --cores=16
```

Prepared input protocols record all row counts, coordinate labels, native levels,
seeds, checksums, calibration indicators and fixed test-U availability. The
preparation tool uses C/X-only folds and proportional native-C stratification
(default C1) with a seeded within-stratum order. O20 is nested in O50. Neither Y
nor physical U participates in choosing folds/masks. Existing supplied folds
are preserved. This is not an unstratified random calibration sample.

The R preparation tool accepts an already quality-controlled, matched cohort;
it does not download raw ADNI records, infer visit/QC rules or embed personal
filesystem paths. It refuses to overwrite prepared inputs.

## Reporting

The seven report sections match Ocean: computation; data/observation design;
numerical assessment and scientific targets; response prediction; response
functions; training-input imputation; reproducibility. ADNI has a multivariate U
and native ordinal ratings, so physical scales and support checks are adapted
explicitly rather than importing Ocean's imposed oxygen bins.

The C definition table lists every model column, ECog field, item meaning and
native grades in the supplied order. Optional `C_labels` in dataset metadata
can describe additional items. Observed-data figures show Y/U and Y/X relations
and the first two C columns against every U coordinate; illustration items are
fixed by input order, not chosen by association strength. All declared C items
remain model inputs. Native grade meanings and unknown-response handling are
stated beside the table.

Tables retain RMSE, normalized RMSE/MAE/CRPS, 95% coverage, normalized interval
width and interval score, fold counts, CORE and native optimizer status. OI and
matched CC are separate. MI is retained in raw outputs but omitted from display.
Unconverged methods remain eligible for numerical bolding and carry an asterisk;
unexposed optimizer status is not called converged.
HTML tables reuse the public Ocean renderer: shaded headers, padded cells,
row rules and horizontal overflow. Entries use mean (sample SD); only the
best mean is bold. The page has no author byline. The table beneath the target
traceplot summarizes that representative fold; full-fold target diagnostics and
integration-node checks are supplementary. MCSE measures numerical uncertainty
in a posterior mean and is distinct from prediction error and posterior SD.

Reports include fixed representative and all-fold U plots, common-hidden scores,
response-free prospective U checks, conditional m/f slices, raw/target traces,
Rhat/ESS/MCSE and latent-integration checks. Derived caches are keyed by fit,
input, numerical settings and analysis source hashes. No report uses held-out Y
to choose an illustrative fold or functional reference profile. Functionals use
64 evenly spaced retained states per chain and256 latent nodes by default; the
fixed target is repeated at512 nodes. These report approximations are separate
from fitting and from full-chain convergence checks. Use
`--functional-draws-per-chain` and `--n-latent` to change them in a new report root.

A compatible saved run is the default report source. A separate local
`--report-source=FILE` descriptor can name historical input/output roots,
a final saved stage and individually versioned comparator recovery locations.
It does not rewrite old fitting identities or refit historical models. The
committed `adni.html` is a data-free plan page; empirical HTML and any individually
named local report variants must not be committed.

The local `reports/paper/` directory also contains ten PDF figures, seven
primary tables in CSV/HTML, and a checksum/status manifest. Figures and tables
are generated from the same saved report analysis.

## Environment

Install the current checkout using the Ocean instructions. The loader checks
installed module hashes against this checkout, rather than trusting a version
number alone. Use `EIVGP_R_LIBRARY` for a custom package library. Dependencies
include jsonlite, digest, ggplot2, patchwork, posterior, knitr, rmarkdown, kergp,
LVGP and EzGP. HTML requires Pandoc on PATH or `RSTUDIO_PANDOC`.

macOS/Linux use the public fork/shared-budget implementation. Windows follows
the package's serial fallback. Use the same input and output flags on all three
entrypoints. Logs, receipts and source provenance stay beside the local results;
failures and nonconvergence are disclosed, not silently reclassified.


## Optional scientific-target ESS and full-history traces

`--ess-min=100` on `adni-eivGP.R` enables saved-state continuation after the
initial sampling budget until every monitored scientific bulk and tail ESS is
strictly greater than 100. The default is `NULL`; `--ess-min=NULL` keeps the
fixed budget. Four chains are required, so this option is unavailable for
runs with fewer than four allocated cores. No new warmup or thinning is introduced. Checks use a
fixed five-point test-input panel (without test responses), 64 fixed latent
integration samples, and all retained MCMC draws. Calibrated cells additionally
monitor physical-surface targets. ESS is a precision criterion, not a proof of
convergence; R-hat and the traceplots remain visible.

The final fit retains `scientific_ess` diagnostics and history. Atomic temporary
continuation states permit resuming an interrupted extension; the final fit is
saved before redundant temporary states are removed. Reporting alone never
starts continuation. With four-chain fits, full-history scientific traces are
also saved when `ess_min` is NULL, without changing the fixed sampling budget.
Two figures show up to four weak mean/uncertainty targets. Display subsampling
spans the complete retained history; reported ESS uses every draw.

ADNI applies the rule only after the final planned stage. Stage directory names
retain the nominal requested budget; convergence tables record actual per-fit
iterations, and the report permits different final iteration counts with the
same warmup and chain count. Figures appear in the saved-results Rmd/HTML.

Install source changes only after active runs have finished. Do not change a
running process or overwrite its installed library. New runs record the ESS
setting in their workflow identity; incompatible historical roots are not
silently relabeled.

For `--manifest=...`, declare `ess_min` in the manifest; any CLI `--ess-min`
value must match it. This avoids silently changing a frozen workflow identity.


### Memory and CPU scheduling

Specify only `--cores=N`. New eivGP fits use `min(N,4)` chains and
`floor(N / min(N,4))` concurrent jobs, capped by pending work. Thus 2 cores
run one two-chain fit, 6 cores run one four-chain fit, 8 cores run two, and
12 cores run three. Prefer multiples of four above four cores. The same
resident-fit limit applies during evaluation. Competitors use one core per
job, allowing up to N concurrent jobs. Ocean follows the same policy.

On a 24 GB Mac, start with `--cores=4`; reduce to 2 if necessary, which changes
the fit to two chains. No automatic RAM guarantee is provided. Changing the
chain count requires a new output directory; saved reports retain the original
chain count. Scientific ESS stopping requires four chains. The former
`--fit-jobs` and `--evaluation-jobs` options are retired. Continuation still
loads each checkpoint once and validates it before use.

Application provenance hashes the ADNI scripts separately from package modules.
Existing results created before this refactor retain their original identity;
source identity checks are not bypassed. Use a fresh output directory for a new
run; migrating existing results requires a separate compatibility audit.
