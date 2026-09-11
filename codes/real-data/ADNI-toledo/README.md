# ADNI comparative case study

`ADNI_Toledo_EIVGP.Rmd` is the analysis entry point. It retains the collaborator's
frozen data, repeats 2/3, four-chain sampler, and checkpoint schedule. Added
comparison/reporting code is in `adni_case_study_helpers.R`. The installed
**eivGP 0.3.1** package performs all EIV fitting, prediction, and imputation.

## Mac mini

See [../MAC_MINI.md](../MAC_MINI.md). From the repository root:

```sh
bash codes/real-data/run_mac.sh check
bash codes/real-data/run_mac.sh smoke --cores 2
bash codes/real-data/run_mac.sh adni --cores 4
```

The launcher uses two workers for four chains and keeps the Mac awake. Set
`--cores 4` to use up to four chain workers. The same flag is supported by
`run_application.R` on Linux. Production runs resume compatible
checkpoints. The `adni` action explicitly disables smoke mode; `smoke` enables it.

## Linux commands

From the repository root, after obtaining the collaborator's private cleaned data:

```sh
Rscript codes/real-data/setup.R   # first setup: also installs kergp, LVGP, EzGP
export ADNI_DATA_DIR="$PWD/real-data/adni"
export ADNI_OUTPUT_DIR="$PWD/results/adni"
EIVGP_SMOKE_TEST=1 EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni
ADNI_REPEAT_ID=2 EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni
# Run repeat 3 on another machine or after repeat 2:
ADNI_REPEAT_ID=3 EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni
```

The short smoke run checks paths and plotting with explicitly reduced competitor
budgets; optimizer failures can occur. It does not establish model performance.
Production competitors use their installed public adapters' default controls.
All packages must be installed before fitting. Code fingerprints and package
versions are saved because 0.3.1 development builds may differ internally.

After transferring both complete repeat output directories into ADNI_OUTPUT_DIR:

```sh
Rscript codes/real-data/ADNI-toledo/report_adni.R
```

This regenerates reports without fitting and writes `combined/table1_prediction.md`
and CSV. It requires complete common-fold comparisons in repeats 2 and 3. Repeat
metrics are averaged through participant losses, with no independent-replicate SE.
The per-repeat figure uses its prespecified fold 1; do not select a favorable
repeat or fold for the paper. Use repeat 2 as the main illustration and repeat 3
as the appendix stability check unless the protocol is amended before results.

## Compact main text

1. Background: PET amyloid burden and fluid biomarkers; explain that ordinal
   proxies are constructed and CSF is naturally unavailable for many participants.
2. `main/fig1_exploration.pdf`: CSF availability, observed-CSF/PET relationship,
   and outcome variation within joint ordinal patterns.
3. `main/table1_prediction.md`: common out-of-fold response prediction for EIV-GP,
   UC-GP, LVGP, and EzGP. It reports RMSE, CRPS, 95% coverage/width, and interval
   score, overall and in the naturally missing-CSF subgroup.
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
