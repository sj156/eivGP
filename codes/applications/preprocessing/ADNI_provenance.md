
## Purpose and data-use boundary

This notebook reconstructs the participant-level inputs used by the
Toledo-inspired ADNI EIV-GP analysis. It creates, **only in a local secure
directory**:

1. `toledo_adni_cohort_n495.csv`;
2. `toledo_adni_balanced_repeated_3fold.csv`;
3. aggregate validation and provenance files.

The generated CSV files contain ADNI participant-level or derived
participant-level data. Each investigator who
runs this notebook must have approved ADNI access and must comply with the
[ADNI Data Use Agreement](https://adni.loni.usc.edu/wp-content/themes/adni_2023/documents/ADNI_Data_Use_Agreement.pdf).


## Raw ADNI files required

Download the following tables through the investigator's own approved
[LONI IDA/ADNI account](https://adni.loni.usc.edu/data-samples/adni-data/).
The displayed download filename may have a date suffix. The notebook searches
recursively by the stable table name, so a file such as
`UCBERKELEY_AMY_6MM_01Aug2026.csv` is accepted.

| Stable table name | ADNI area / description | Columns used |
|---|---|---|
| `UPENN_PLASMA_FUJIREBIO_QUANTERIX` | Biofluid Biomarkers; UPENN plasma Aβ42, Aβ40, p-tau217, NfL and GFAP measured by Fujirebio/Quanterix | `RID`, `EXAMDATE`, `pT217_F`, `pT217_AB42_F`, `AB42_AB40_F`, `GFAP_Q` |
| `UCBERKELEY_AMY_6MM` | PET numerical data; UC Berkeley amyloid PET 6 mm analysis | `RID`, `SCANDATE`, `CENTILOIDS`, `SUMMARY_SUVR`, `TRACER`, `qc_flag` |
| `DXSUM` | Clinical Assessments; diagnostic summary | `RID`, `EXAMDATE`, `DIAGNOSIS` |
| `PTDEMOG` | Subject Characteristics; participant demographics | `RID`, `VISDATE`, `PTGENDER`, `PTEDUCAT`, `PTDOB` |
| `APOERES` | Genetics / Biospecimen; ApoE genotyping results | `RID`, `GENOTYPE` |
| `UPENNBIOMK_ROCHE_ELECSYS` | Biofluid Biomarkers; UPENN CSF Roche Elecsys biomarkers | `RID`, `EXAMDATE`, `ABETA42`, `BATCH` |

For the exact frozen experiment, use the same ADNI snapshot used in the
original analysis (the PET, diagnosis, and demographics exports were dated
01-Aug-2026). ADNI tables are updated over time. A newer snapshot can change
the eligible cohort and therefore cannot silently replace the frozen
confirmatory cohort.

The preferred input is the original CSV exported from IDA. For local backward
compatibility, this notebook can also read a matching `.rda` file or a matching
CSV stored inside a ZIP archive. Those compatibility options do not authorize
redistribution of the data.

## Configuration

Before rendering, set two environment variables. Example:

```r
Sys.setenv(
  ADNI_RAW_DIR = "/secure/path/to/my/adni-downloads",
  ADNI_TOLEDO_DATA_DIR = "/secure/path/to/my/adni-derived/toledo"
)
rmarkdown::render("ADNI_preprocess.Rmd")
```

Raw downloads default to `real-data/raw/adni/`; cleaned output defaults to
`real-data/` at the repository root. This private directory is excluded from Git.
The historical repeated-fold file produced below is retained for provenance;
`ADNI_cv.R` generates the current five-fold analysis assignment.

Set `ADNI_STRICT_FROZEN=0` only for an explicitly labelled cohort-update
sensitivity analysis. Such output is not the frozen n=495 experiment.


## Input discovery and validation helpers

The environment-variable overrides below are optional. They are useful if the
raw directory contains multiple versions of the same ADNI table:

```text
ADNI_PLASMA_FILE
ADNI_PET_FILE
ADNI_DIAGNOSIS_FILE
ADNI_DEMOGRAPHICS_FILE
ADNI_APOE_FILE
ADNI_CSF_FILE
```

Each override should point to the exact local CSV, RDA, or ZIP file. When a ZIP
is supplied, the notebook selects the member beginning with the corresponding
stable table name.


## Load the six authorized raw tables

No participant rows are printed in this report.


## Build the Toledo analysis cohort

The construction order reproduces the original frozen analysis:

1. use the plasma visit as the temporal anchor;
2. match the nearest QC-passed amyloid PET scan within +/-183 days;
3. match the nearest diagnosis to PET within +/-183 days;
4. attach demographics and APOE4 dose;
5. retain the closest eligible plasma-PET landmark per RID **before** reading
   CSF availability;
6. define observed `U` as Roche Elecsys CSF Aβ42 within +/-30 days of PET;
7. retain florbetapir (`FBP`) PET only;
8. derive plasma Aβ42/Aβ40 and GFAP tertiles from the complete n=495 cohort
   without using `Y`, `U`, or `R`.


## Construct the frozen balanced 3 repeats x 3 folds

This is the exact balance-only procedure used for the current `balanced-v2`
split. It never fits a prediction model and never selects a split based on
held-out predictive performance.

Each fold has exactly 165 participants. Balance targets include diagnosis,
natural U availability, outcome range, sex, APOE4 dose, the two-proxy joint
cells, age, observed-U range, and coarse observed-U/amyloid-PET curve regions.
The use of `Y` and observed `U` here is only for allocation balance; it does not
change cohort membership and does not evaluate a fitted model.


## Final frozen-reproduction check and local manifest

For the 01-Aug-2026 source snapshot, the expected output hashes are included as
reproducibility checks. A hash mismatch is not repaired by editing the cohort or
fold files manually. First check raw-table versions, column definitions, locale,
and R version.


## Session information

The report records the local software environment but never prints participant
rows.

