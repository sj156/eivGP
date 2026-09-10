# ADNI application and ocean real-data validation

Both examples use **eivGP exactly 0.3.1**, through the public package fitting
and prediction APIs. The notebooks are the runnable analysis entry points:

- `ADNI-toledo/ADNI_Toledo_EIVGP.Rmd`: observational ADNI application;
  multivariate ordinal-probit model, naturally missing CSF input.
- `Ocean_Validation.Rmd`: validation with controlled hiding of measured
  phosphate; univariate deterministic threshold model.

## Linux: pull, then run

From the eivGP repository root:

```sh
git pull
# First machine setup only: installs dependencies and tested eivGP 0.3.1 commit.
Rscript codes/real-data/setup.R

# Short ADNI pipeline check:
EIVGP_SMOKE_TEST=1 EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni
# Production ADNI (repeat 2; folds 1–3):
ADNI_REPEAT_ID=2 EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni
# Repeat 3 can run on another computer:
ADNI_REPEAT_ID=3 EIVGP_N_CORES=4 Rscript codes/real-data/run_application.R adni

# Ocean requires its prepared data:
OCEAN_DATA_DIR=/path/to/prepared OCEAN_QUICK=true MIXEDGP_CORES=4 Rscript codes/real-data/run_application.R ocean
# Replace OCEAN_QUICK=true with false for production.
```

`run_application.R` extracts and executes the selected Rmd's R chunks. It
requires knitr but neither RStudio nor Pandoc. For an HTML report instead,
install rmarkdown and Pandoc and run:

```sh
Rscript -e 'rmarkdown::render("codes/real-data/ADNI-toledo/ADNI_Toledo_EIVGP.Rmd")'
Rscript -e 'rmarkdown::render("codes/real-data/Ocean_Validation.Rmd")'
```

Setup pins commit `a3aa0f697f97d4bdf027507ac3c9b952644affb9`. Runtime requires
version 0.3.1 and accepts local builds; ADNI records the installed commit when
available. Transitive dependency versions are not locked. No package is
installed automatically during analysis.

## Data size and distribution

ADNI cohort: **495 rows × 29 columns**, 75,441 bytes as CSV. There are 133
observed and 362 naturally missing CSF inputs. The repeated-fold file has
**1,485 rows × 3 columns**, 13,612 bytes (three assignments per participant).
Each fold has 330 training and 165 test participants. Combined CSV size is
89,053 bytes (about 87 KiB).

Ocean's allocation table specifies **100 training observations** across six
classes. The prepared `study1_data.csv` and `meta.json` are absent from this
repository; total rows, test rows, and dataset bytes remain unverified.

R package data would normally load with `data("ocean", package="eivGP")`,
not `load()`, which loads a file by path. Such a dataset is not yet implemented.
Ocean requires the source files and verification of redistribution terms.
The ADNI DUA prohibits redistribution of participant-level data, including
derived data, so the actual cohort should not be shipped as public package
data. See https://adni.loni.usc.edu/wp-content/themes/adni_2023/documents/ADNI_Data_Use_Agreement.pdf.
An independently generated synthetic example is an option for a public ADNI-style
package demonstration, clearly distinguished from the actual study.

Authorized ADNI users can set `ADNI_DATA_DIR=/private/path` containing the two
frozen CSVs. The fallback remains `ADNI-toledo/data/` for compatibility with the
collaborator's existing local folder. Data hashes must match the frozen design.
The CSVs already tracked in this repository have not been removed by this
cleanup; their presence does not establish permission for public distribution.

## Outputs and resuming

ADNI outputs default to `ADNI-toledo/outputs/repeat2/` or `repeat3/`; smoke
outputs use `smoke_repeat2/` etc. `ADNI_OUTPUT_DIR` overrides the output base.
Four chains run in parallel; folds run sequentially. Production uses 20,000
transitions (4,000 warm-up), with one extension to 30,000 if diagnostics fail.
Checkpointing begins at iteration 6,000, then occurs every 2,000 transitions.
Rerun the same command to resume. Copy the entire repeat output directory to
move a checkpoint between computers. Configuration mismatches stop the run.
A repeat-level `.run-lock` prevents simultaneous writers. Remove a stale lock
only after confirming its process has stopped. Actual iteration counts label
outputs; `convergence_passed` distinguishes finished execution from a passing
convergence gate. Smoke draws are not suitable for inference.

Ocean outputs default to `outputs/`. `OCEAN_OUTPUT_DIR`, `OCEAN_KERNEL=se|matern`,
`OCEAN_MATERN_NU=2.5`, and `OCEAN_USE_CACHE=true|false` are configurable. Its
cache schema separates the migration to the 0.3.1 API from old caches. Use a
new output directory when data change or for concurrent runs. Ocean has no
within-chain resume. Shared historical baseline and reporting helpers are
bundled under `shared/`; EIV fits and response prediction use the package.

Ocean predicts from X and ordinal C without test U. ADNI uses observed test U
where available and integrates missing U, reporting scores overall and by
missingness. The supplied ADNI notebook fits EIV-GP only and permits repeats
2 and 3, preserving the collaborator's experiment design.

## Raw-data processing and provenance

Read [Raw_Data_Processing.Rmd](Raw_Data_Processing.Rmd) for the known processing
steps, evidence sources, and missing original procedures. Validate prepared
inputs without starting MCMC:

```sh
Rscript codes/real-data/run_application.R data
```

The original raw-to-prepared scripts are still needed from the collaborator.
This command checks existing prepared data; it does not reconstruct them.
The manuscript specifies 100 training plus 300 test ocean records; this updates
the size description above, but the prepared files remain unavailable to verify.
