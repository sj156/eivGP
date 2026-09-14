# Real-data applications

Run from the repository root on macOS or Linux. Both applications require
**eivGP 0.3.1**. The R runner executes the analysis notebooks with knitr; RStudio,
Pandoc, and shell launchers are not required.

```text
eivGP/
  codes/
    ADNI_application.Rmd          # Scientific motivation, fitting and results
    ADNI_preprocess.Rmd           # Collaborator raw-to-cleaned cohort workflow
    ADNI_cv.R                    # Fixed five-fold assignment
    ADNI_report.R                # Rebuild the pooled prediction report
    ADNI_case_study_helpers.R     # Competitors and reporting
    ADNI_parallel.R              # Core-budget scheduler
    ocean_validation.Rmd         # Controlled-hiding validation
    ocean_preprocess.Rmd         # Collaborator raw-to-cleaned ocean workflow
    ocean_01_representative_figures.R
    ocean_02_replicates.R
    ocean_run_all.R
    ocean_data_helpers.R
    run_real_data.R               # Common Rscript entry point
    real_data_paths.R             # Data and output defaults
    real_data_checks.Rmd          # Prepared-data checks, no fitting
    setup_real_data.R
    check_real_data.R
    real_data_support/            # Shared ocean support code
    tests/ADNI_*.R
  eivGP/                         # R package source; applications are separate
  real-data/                     # Private prepared data, excluded from Git
    toledo_adni_cohort_n495.csv
    ocean/
      study1_data.csv
      meta.json
      eligible_pool.csv
  real_data_application_outputs/ # Generated files, excluded from Git
    adni/five_fold/               # main/, appendix/, fold_1/ ... fold_5/
    adni/smoke/                   # Short pipeline-check results
    ocean/                       # figures/, tables/, fits/
```

The previous `codes/real-data/` tree has been flattened into `codes/`.
Code alone can be pulled with Git; private inputs must be transferred separately
by an authorized user. Results may contain participant-level information and
are also excluded from Git. Absolute data/output overrides remain available:
`ADNI_DATA_DIR`, `OCEAN_DATA_DIR`, `ADNI_OUTPUT_DIR`, `OCEAN_OUTPUT_DIR`.

## Run ADNI

```sh
# First machine setup only (requires internet):
Rscript codes/setup_real_data.R
Rscript codes/check_real_data.R

# Allocation preview: no data loaded or fitting started.
Rscript codes/run_real_data.R adni --cores 12 --plan

# Optional short pipeline check, isolated from production:
EIVGP_SMOKE_TEST=1 Rscript codes/run_real_data.R adni --cores 12

# Production: all five folds and the pooled comparison report.
EIVGP_SMOKE_TEST=0 Rscript codes/run_real_data.R adni --cores 12

# Rebuild reports later, without fitting:
Rscript codes/ADNI_report.R
```

No repeat or validation arguments are needed. Old `--repeats` and `--validations`
options fail with a migration message. Old environment settings are ignored
with a warning. The historical three-fold CSV is no longer an input dependency.

The fixed seed is 20260914, using explicit R RNG settings. Participants are
sorted by RID and shuffled within observed/missing-CSF strata. Each test fold
has 99 people, including 26 or 27 with measured CSF; training uses the other
396. Neither PET outcomes, measured CSF values nor method performance choose
the split. `appendix/fold_assignments.csv` records the exact assignment.

Five-fold CV trains on four groups and rotates the held-out group through all
five. This gives one held-out prediction per participant, with five fits per
method. Ten-fold CV would use 445–446 training and 49–50 test participants and
require ten fits per method. Five folds is a practical choice here because GP
MCMC is costly. Although there are fewer fits than the old three-times-three
scheme, each uses more training data; a shorter total runtime is not guaranteed.
Do not change or choose splits based on the resulting method rankings.

The main table pools participant losses across all 495 held-out predictions;
RMSE is the square root of the pooled mean squared error. It does not average
fold RMSEs or treat folds as independent studies. Main figures and the comparison
table are in `main/`; subgroup metrics, failures, paired losses, hidden-CSF
validation, package versions and diagnostics are in `appendix/`. The physical
CSF surface uses prespecified fold 1. All methods use identical train/test sets;
EIV-GP's additional training CSF information is explicitly reported.

## Resources, checkpoints and direct notebooks

`--cores N` is the total worker budget. All fits retain four MCMC chains.
Examples: 12 cores -> three concurrent folds x four chain workers; 16 -> four x
four; 56 -> five x four (20 active chain workers maximum for this design).
Unused capacity does not create extra folds or chains. Comparator phases and
reporting can use fewer cores. Set `ADNI_MAX_FOLD_WORKERS` to limit memory use.

Production keeps 20,000 iterations per chain, 4,000 warm-up, and one extension
to 30,000 when diagnostics fail. A completed run is not necessarily converged;
check diagnostics. Repeat the same command to resume compatible checkpoints.
The new five-fold design cannot resume old three-fold or numbered-validation
fits. A run lock and per-fold locks prevent concurrent writers. Remove a stale
lock only after its process has stopped. Inspect `fold_N/worker.log` and
`CURRENT_STATUS.txt` while running. No analysis starts when using `--plan`.

To render the notebook itself (rmarkdown and Pandoc required), write its HTML
outside codes as well:

```sh
Rscript -e 'd <- "real_data_application_outputs/adni"; dir.create(d, recursive=TRUE, showWarnings=FALSE); Sys.setenv(EIVGP_CORES=12); rmarkdown::render("codes/ADNI_application.Rmd", output_dir=normalizePath(d))'
```

This runs the same five folds. Do not run it alongside the Rscript coordinator.
For macOS sleep prevention, prefix the production command with `caffeinate -i`.
The optional `run_real_data_mac.sh` wrapper remains available; it is not required.
See `REAL_DATA_MAC.md` for the compact Mac commands.

## Run ocean

```sh
# Check both prepared datasets without fitting:
Rscript codes/run_real_data.R data
# Pipeline check; production uses OCEAN_QUICK=false:
OCEAN_QUICK=true Rscript codes/run_real_data.R ocean --cores 12
```

Ocean retains its fixed 100 training / 300 test design and nested controlled
hiding of training phosphate; ADNI's five-fold choice does not alter it.
The supplied metadata records calibration counts per class. The loader uses
those if the optional `class6_proportional_calibration_allocation.csv` is absent;
it never estimates them from outcomes. Ocean fits use cache reuse rather than
within-chain resume. Use a separate output override when changing data or
settings. Set `OCEAN_USE_CACHE=false` for an intentional recomputation.

## Preprocessing and provenance

`ADNI_preprocess.Rmd` and `ocean_preprocess.Rmd` are the collaborator's raw-data
workflows. Their cleaned outputs default to `real-data/` and `real-data/ocean/`.
ADNI raw downloads default to `real-data/raw/adni/`; ocean defaults to
`real-data/raw/ocean/AllBottle.csv`. Raw-file access and cohort construction
remain separate from fitting. `real_data_checks.Rmd` validates prepared inputs;
it does not reproduce the original raw-data processing.

The frozen ADNI cohort has 495 rows and 29 columns, 133 observed and 362 missing
CSF measurements. Its checksum is recorded in `ADNI_DATA_MANIFEST.txt`.
Upstream cohort selection, assay coarsening and previous development decisions
still require provenance review; a fresh CV assignment does not undo leakage
or model selection. Sampler settings, package version requirement, scientific
tasks and the published competitor set are unchanged by this structural cleanup.
