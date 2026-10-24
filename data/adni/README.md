# Local ADNI data

This folder contains a data contract and metadata example, **no ADNI participant
records**. Authorized collaborators supply the prepared files locally after
cloning. CSVs and private JSON protocols are ignored by Git; protocols include
participant IDs and observation masks.

## 1. Supply verified prepared inputs

Place the files directly in `data/adni/`, without `inputs/`, `MC/` or per-fold
subdirectories:

```text
data/adni/
  DATASET.json
  fold01_O0_PROTOCOL.json
  fold01_O0_train.csv
  fold01_O0_test.csv
  fold01_O20_PROTOCOL.json
  fold01_O20_train.csv
  fold01_O20_test.csv
  ...                         # through fold05_O50
```

There are 45 fold files (five folds × three calibration cells × three files)
plus `DATASET.json`. The separately supplied three-item cohort (MEMORY3, MEMORY8,
DIVATT3) uses the same filenames directly in the sibling folder `data/adni_3C/`.
Add `--dataset=3C` to all three application commands; paper results go to
`results/3C/`. The default nine-item cohort uses `data/adni/` and `results/adni/`.
Both input folders are excluded from Git except this contract and the example.

Each CSV contains `RID,Y,X1,...,U1,...,C1,...,calibration`. Numeric X columns are
continuous. U reference values are retained for scoring; the fitter masks every
noncalibration training U vector. All test calibration values are zero. Native
ordinal C items retain integer grades `1:m_j`. Participants are held out once
across five disjoint test folds. O0/O20/O50 share the same rows and values;
O20 calibration and test-U availability are nested in O50.

Protocols supply the dataset identifier `plan_id`, `p,q,d`, coordinate labels,
`m_vec`, `C_levels`, fold/cell/counts, `ident`, `standardize_U`, fixed
`test_U_available_RID`, CSV MD5s and task seeds. O0 uses
`ident="lower_triangular",standardize_U=false`; O20/O50 use
`ident="none",standardize_U=true`. Physical calibration vectors must have
sufficient affine rank. Keep the supplied CSV/JSON bytes unchanged.

Older prepared trees (`MC/foldNN/O*/private/*.csv` and `PROTOCOL.json`) remain
readable through `--data-root`; new preparation and result snapshots use flat
filenames. Existing frozen data are not regenerated or moved automatically.

## 2. Alternatively, prepare a local matched cohort

Create a private CSV with one eligible row per participant, the columns above
without `calibration`, and optionally a previously frozen `fold` column. Create a
private JSON using `DATASET.example.json`, adapting coordinate labels, native
level counts and dataset identity. Different C sets require matching `C_order`,
`m_vec` and `C1,...,Cq` columns; no R code changes are needed.

```sh
Rscript tools/adni-prepare-folds.R /private/cohort.csv /private/DATASET.json data/adni
```

The caller supplies an eligible, matched, quality-controlled cohort. Without
supplied folds, preparation balances C levels and continuous X only.
Calibration/test-U sampling uses fixed seeds within native-C strata (C1 by
default) with proportional nested 20%/50% quotas. Neither Y nor U chooses folds
or masks. Preparation refuses to overwrite frozen fold files. Changing seeds or
the preparation algorithm creates a different frozen partition.

## 3. Run and keep outputs local

Use the three commands in [the application guide](../../paper/applications/README_adni.md).
`--data-root=/external/adni` reads files outside the checkout. New result roots
keep a verified flat copy in `results/adni/data/`. Fits and empirical reports
remain ignored. Share participant inputs and derived reports through an
appropriate authorized channel; the repository distributes portable code,
the contract and a data-free HTML plan.
