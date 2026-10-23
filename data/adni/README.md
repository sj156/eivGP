# Local ADNI data

This folder contains code and a metadata example, **no ADNI participant records**.
Git cannot retain an empty directory, so `inputs/.gitkeep` is only a placeholder.
Authorized collaborators must supply the data locally after cloning. All files
beneath `inputs/` except that placeholder are ignored by Git, including JSON
protocols, which contain participant IDs and observation masks.

## Option A: copy a verified prepared input directory

Place the following 45 files under `data/adni/inputs/` (five folds × three cells
× three files per cell):

```
inputs/MC/fold01/O0/PROTOCOL.json
inputs/MC/fold01/O0/private/train.csv
inputs/MC/fold01/O0/private/test.csv
... fold01/O20, fold01/O50, ... fold05/O50
```

Each CSV contains `RID,Y,X1,...,U1,...,C1,...,calibration`. Numeric X columns are
continuous. U reference values are retained for scoring, while the fitter masks
every noncalibration training U vector before calling the model. All test
`calibration` values are zero. C items use native integer grades `1:m_j`; no
new imaging/response bins are constructed. Participants are held out once across
five disjoint test folds. O0/O20/O50 share the same rows and values; only training
calibration differs. O20 training/test masks are nested in O50.

Protocols supply `plan_id` (a portable dataset identifier), `p,q,d`, `X_order`,
`U_order`, `C_order`, `m_vec`, `C_levels`, fold/cell/counts, `ident`,
`standardize_U`, fixed `test_U_available_RID`, CSV MD5s and task seeds.
O0 uses `ident="lower_triangular",standardize_U=false`; O20/O50 use
`ident="none",standardize_U=true`. All physical calibration vectors must have
sufficient affine rank. Use `prepare_folds.R` to generate the complete contract.

The separately supplied three-item cohort (MEMORY3, MEMORY8 and DIVATT3) uses
`data/adni/3C/inputs/` with the same structure. Run all three application
entrypoints with `--dataset=3C --profile=paper`; outputs go to
`results/3C_paper/`. All participant inputs in this additional folder are also
ignored by Git.

## Option B: prepare from a local matched cohort

Create a private CSV with one eligible row per participant, the columns above
without `calibration`, and optionally a previously frozen `fold` column. Create a
private JSON using `DATASET.example.json`, adapting coordinate labels, native
level counts and dataset identity. The example is the nine-item memory/attention
setting. If using a different number of C items, provide matching `C_order` and
`m_vec` arrays and `C1,...,Cq` columns; no R code changes are needed.

```sh
Rscript data/adni/prepare_folds.R /private/cohort.csv /private/DATASET.json data/adni/inputs
```

The caller is responsible for ADNI permissions, female restriction, measurement
matching, native grades and source QC before preparation. The tool does not
select or filter participants by Y/U performance. Without supplied folds it
balances C levels and continuous X only. Calibration/test-U sampling uses fixed
seeds within native-C strata (C1 by default) and proportional nested20%/50%
quotas. It never overwrites an existing `MC` directory. A different implementation
or seed is a new frozen partition, not the same historical numerical realization.

## Run and keep outputs local

Use the three commands in [the application guide](../../paper/applications/README_adni.md).
`--data-root=/external/inputs` avoids copying data into the checkout at all.
`results/` and empirical HTML remain ignored. Share participant inputs and derived
reports only through an appropriate authorized channel; the repository contains
only the portable code/contract and a data-free HTML plan.
