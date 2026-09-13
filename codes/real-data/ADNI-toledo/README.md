# ADNI comparative case study

`ADNI_Toledo_EIVGP.Rmd` is the analysis entry point. It retains the collaborator's
frozen data, four-chain sampler, and checkpoint schedule. Added
comparison/reporting code is in `adni_case_study_helpers.R`. The installed
**eivGP 0.3.1** package performs all EIV fitting, prediction, and imputation.

## Twenty-validation split extension

The 20-validation protocol does not duplicate the private cohort into twenty
analysis datasets. It keeps one cohort and creates twenty small assignment
files containing only `RID`, `train`/`test` membership, and reproducibility
metadata. Generate them locally with:

```sh
export ADNI_DATA_DIR="$PWD/real-data/adni"
Rscript codes/real-data/ADNI-toledo/prepare_validation_splits.R
```

The files are written to `ADNI_DATA_DIR/validation-splits/` unless
`ADNI_SPLIT_DIR` is set. `validation_01.csv`--`validation_09.csv` copy the
existing frozen 3 repeats x 3 folds exactly. The original three repeat-level
assignment seeds are recorded truthfully: the three folds within a repeat were
generated jointly and therefore do not have independent split seeds.
`validation_10.csv`--`validation_20.csv` use the frozen seeds in
`validation_split_plan.csv` to draw separately randomized one-third test holdouts within
observed-CSF status. The only rejection rule preserves at least two
observed-CSF training participants in every joint proxy cell; outcome values
and fitted-model performance are not used to select a new split.

The generated files contain participant identifiers and must be transferred
only through approved private storage, never committed to this public
repository. The validations reuse and overlap participants, so the reported SD
is a descriptive stability summary.

### Preferred 20-validation workflow

The public repository contains the generator and frozen seed registry, not the
private assignment CSVs. On one trusted machine, generate the assignments once
and privately copy the complete `validation-splits/` directory to every worker
machine together with the two existing cleaned input CSVs. Do not regenerate
or hand-edit individual assignment files on different machines: the private
manifest is checked before fitting.

Each machine can then claim disjoint validation IDs. For example:

```sh
# Machine A: inspect its allocation, then run validations 1--5.
Rscript codes/real-data/run_application.R adni --cores 16 --validations 1-5 --plan
Rscript codes/real-data/run_application.R adni --cores 16 --validations 1-5

# Machine B (same code, cohort, and validation-splits directory):
Rscript codes/real-data/run_application.R adni --cores 16 --validations 6-10

# Individual/non-contiguous IDs are also accepted.
Rscript codes/real-data/run_application.R adni --cores 8 --validations 11,14-16
```

One validation is one 330-person training fit plus one 165-person test set.
Its output is self-contained in `ADNI_OUTPUT_DIR/validation_XX/`; copying that
whole directory does not overwrite any other validation. After collecting all
20 directories under the same `ADNI_OUTPUT_DIR`, rebuild the combined report:

```sh
Rscript codes/real-data/ADNI-toledo/report_adni.R --validations 1-20
```

The combined report writes `combined/all_validation_main_metrics.csv` and
`combined/table1_prediction_20validations.{csv,md}`. Across-validation SDs are
descriptive stability summaries only because test sets overlap.

## Mac mini

See [../MAC_MINI.md](../MAC_MINI.md). From the repository root:

```sh
bash codes/real-data/run_mac.sh check
bash codes/real-data/run_mac.sh smoke --cores 2 --validations 10
bash codes/real-data/run_mac.sh adni --cores 4 --validations 10-12
bash codes/real-data/run_mac.sh report --validations 1-20
```

The launcher defaults to a two-core budget and keeps the Mac awake. Set
`--cores N` to divide a larger budget across validation jobs and chains. The same flag is supported by
`run_application.R` on Linux. Production runs resume compatible
checkpoints. The `adni` action explicitly disables smoke mode; `smoke` enables it.

## Linux commands

From the repository root, after obtaining the collaborator's private cleaned data:

```sh
Rscript codes/real-data/setup.R   # first setup: also installs kergp, LVGP, EzGP
export ADNI_DATA_DIR="$PWD/real-data/adni"
export ADNI_OUTPUT_DIR="$PWD/results/adni"
EIVGP_SMOKE_TEST=1 EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni --validations 1
EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni --validations 1-5
# Run disjoint IDs on another machine:
EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni --validations 6-10
```

The short smoke run checks paths and plotting with explicitly reduced competitor
budgets; optimizer failures can occur. It does not establish model performance.
Production competitors use their installed public adapters' default controls.
All packages must be installed before fitting. Code fingerprints and package
versions are saved because 0.3.1 development builds may differ internally.

For the preferred protocol, collect `validation_01` through `validation_20` in
the same ADNI_OUTPUT_DIR and run:

```sh
Rscript codes/real-data/ADNI-toledo/report_adni.R --validations 1-20
```

This regenerates reports without fitting. It requires complete common-method
comparisons for every requested validation and reports the arithmetic mean and
descriptive SD of each validation-level metric.

## Compact main text

1. Background: PET amyloid burden and fluid biomarkers; explain that ordinal
   proxies are constructed and CSF is naturally unavailable for many participants.
2. `main/fig1_exploration.pdf`: CSF availability, observed-CSF/PET relationship,
   and outcome variation within joint ordinal patterns.
3. `main/table1_prediction.md`: common out-of-fold response prediction for EIV-GP,
   UC-GP, LVGP, and EzGP. It reports RMSE, CRPS, 95% coverage/width, and interval
   score overall and for both R=0 and R=1 subgroups.
4. `main/fig2_CSF_inference.pdf`: held-out masked-CSF validation and a pointwise
   posterior CSF-coordinate response surface. Neither naturally missing CSF nor
   the true response surface is observed, so those cannot be assigned truth-based
   accuracy scores. PNG previews accompany both figures.

The published methods follow the computation section's main method roster and
use the public eivGP adapters to kergp, LVGP, and EzGP. CC-GP and PI-GP remain
separate calibration-aware simulation ablations; they have not been silently
replaced by differently defined ADNI methods. The appendix LM-CSF benchmark is
explicitly a simple response-free imputer, not PI-GP or a published GP method.

## Comparison protocol

- Primary: all methods predict test Y given X,C, with test CSF withheld. EIV-GP
  uses observed training CSF; the literature competitors use X,C,Y. The comparison
  shows utility of the complete workflows, not an isolated algorithmic advantage
  at identical training-information use.
- Secondary: use available test CSF for EIV-GP, retaining X,C predictions for
  categorical competitors. These predictions are reported in the appendix.
- Prospective CSF validation hides observed test CSF, passes only test C to
  `impute_eivgp`, and never passes test Y. It evaluates the observed-CSF subset,
  not naturally missing CSF. The LM-CSF benchmark uses the same training observed
  CSF and test ordinal variables. Training imputation uses training Y and is
  saved separately without unsupported truth scores.
- Surface plot: fold 1, median training age, female=0, APOE4 dose=1, observed
  training-CSF 10th–90th percentile range. It shows pointwise intervals and an
  association, not a causal effect or a population-averaged curve.
- All transformations and competitor tuning use training data. Upstream cutoff
  construction and screening history still need the collaborator's documentation.
- Invalid public-package predictions remain failures. No variance patch, substitute
  implementation, or favorable-fold selection is introduced. Main-table subsets
  are common to all methods and are labelled incomplete unless all folds finish.
- Never interpret smoke results as evidence. Failed MCMC gates and optimizer
  warnings remain visible; do not assert superiority from untrustworthy estimates.

## Appendix and private outputs

`appendix/` contains fold/scenario scores (including MAE, NLPD when conditional
Gaussian components are available, and 50/80/95% coverage), method availability,
paired loss contrasts, CSF validation scores, cohort summaries, and package
fingerprints. Paired bootstrap intervals condition on fitted folds; they are
not full repeated-CV uncertainty intervals and do not include refitting variation.

Each `fold_N/case-study/` contains bounded EIV predictive draws, cached competitor
predictions, participant-level scores, prospective CSF checks, missing-training-CSF
posteriors, all hyperparameter diagnostics, dictionary occupancy, and prediction
scores by chain. Participant-level outputs must remain in private storage.

Reruns reuse the EIV checkpoint and compatible comparator cache. Completed failures
are cached too; inspect them, then remove that fold's `case-study/competitors.rds`
explicitly if a rerun is warranted. Changed comparator settings/input/package
fingerprints reject cache reuse. Full MCMC checkpoint compatibility is unchanged.
`EIV_COMPLETED.txt` marks the original fitting/prediction task;
`CASE_STUDY_COMPLETED.txt` marks completion of the extended workflow. Both still
require inspection of convergence and competitor status.

## Data detection for direct R execution

Both `Rscript codes/real-data/run_application.R adni --cores 4` and direct
rendering of the Rmd detect the two cleaned CSVs in `real-data/adni/`, then
`real-data/`, then the legacy `codes/real-data/ADNI-toledo/data/`. An explicit
`ADNI_DATA_DIR` overrides detection and is checked without falling back to a
different dataset. Missing files now stop before an output lock or MCMC run.

## Parallel validation scheduling

The validation interface schedules separate validation jobs:

```sh
Rscript codes/real-data/run_application.R adni --cores 16 --validations 1-20 --plan
Rscript codes/real-data/run_application.R adni --cores 16 --validations 1-20
```

The coordinator treats every validation as one fit unit and dynamically fills
the available worker slots. For example, 16 cores run four validation fits at a
time with four chains per fit. Follow `validation_XX/fold_1/worker.log`,
`PROGRESS.md`, or `CURRENT_STATUS.txt`. The scheduler status file is named
`scheduler-validation-...-status.csv`.

The coordinator holds each selected validation's `.run-lock`. Each worker holds
its own `.fold-lock` and writes only that validation's files. Folds are retained
as an internal directory label (`fold_1`) so existing checkpoint logic remains
simple. A failure preserves other completed checkpoints, stops combined
reporting, and records its error.

An already-running process retains the code it loaded. Do not start the new
coordinator alongside it. To switch, allow a checkpoint to complete, interrupt
the old run, confirm its R processes/workers have stopped, and rerun with the
new command. Remove a stale validation/fold lock only after confirming its owner is
gone. Existing compatible checkpoints are reused; an interrupted unsaved segment
must be recomputed. Changing concurrency does not change checkpoint compatibility.

Direct Rmd rendering runs one selected validation via `ADNI_VALIDATION_ID`;
`run_application.R` coordinates multiple validation IDs. Smoke outputs remain
separate from production outputs.
