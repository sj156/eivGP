# Local OAI prepared inputs

This folder contains a data contract and a metadata example, no participant
records. Authorized collaborators place privately supplied files here after
cloning. Git ignores participant CSVs, private JSON protocols, raw archives,
model fits and empirical reports.

## 1. Supply a frozen input bundle

Use the same flat format as ADNI:

```text
data/oai/
  DATASET.json
  fold01_O0_PROTOCOL.json
  fold01_O0_train.csv
  fold01_O0_test.csv
  fold01_O20_PROTOCOL.json
  fold01_O20_train.csv
  fold01_O20_test.csv
  ...                          # through fold05_O50
```

There are 45 fold files plus `DATASET.json`. `--data-root=/local/oai` can point to
another directory. For another named cohort, use `data/oai_NAME/` with
`--dataset=NAME`; paper output then goes to `results/NAME/`.

Each CSV contains `RID,Y,X1,...,U1,...,C1,...,calibration`. X is continuous; constant
sex is handled by cohort eligibility rather than a binary GP input. Native
ordinal C codes are 1:m_j, with original-grade mappings supplied in metadata.
Raw U references stay in private files for scoring; fitting masks every U in
noncalibration training rows. Test CSV calibration is zero. Protocols identify
which held-out U may be used for matched CC scoring.

The five test folds are disjoint and share a fixed cohort. O0/O20/O50 have
identical row order and X/C/Y/U values. O20 training calibration and test-U
revelation are nested in O50. O0 uses an unanchored model scale; O20/O50 use
sufficient-rank physical calibration. Protocols include `plan_id`, `p,q,d`,
coordinate orders, `m_vec`, level maps, fold/cell/counts, `ident`, `standardize_U`,
`test_U_available_RID`, CSV MD5s and fixed task seeds. Keep supplied bytes fixed.
Older `MC/foldNN/O*/private/` trees remain readable, but new preparation and
result snapshots are flat.

The default metadata describes baseline V00JSW225, baseline V00WMTFMTH, age/BMI/
height and eight native osteophyte/sclerosis grades. Quantitative MRI, quantitative
radiographic JSW and semiquantitative ordinal readings come from distinct OAI
assessment projects. No data are downloaded by the application.

## 2. Prepare a new eligible cohort

Provide a private cohort CSV with the same variables except `calibration`, plus
an optional frozen `fold` column. Adapt `DATASET.example.json` for the actual
native C mappings, coordinate definitions and fixed seeds:

```sh
Rscript tools/oai-prepare-folds.R /local/cohort.csv /local/DATASET.json data/oai
```

Preparation uses C/X-only balanced folds and proportional
random within-stratum nested 20/50 masks. The default OAI example uses joint
C1/C2 strata. Y/U do not choose folds or masks. Existing prepared folds are never
overwritten; no rare level is silently merged. Input validation requires each
native grade in every full training fold. Check calibration-category support
separately for small CC training sets.

A newly prepared cohort is a new allocation. For reproduction, use the supplied
frozen bundle unchanged rather than running preparation again. All methods read
the same training-calibration and test-U masks directly from its protocols.

## 3. Keep execution records local

Follow [the OAI guide](../../paper/applications/README_oai.md). The default paper
root is `results/oai/`; it retains verified input snapshots and saved reports.
Participant inputs, masks, raw downloads and derived fits are not committed.
