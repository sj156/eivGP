
This notebook is the reproducible data-preparation record for the Ocean
application in the EIV-GP experiments. It reads the public BCO-DMO
Dataset 3228 file `AllBottle.csv`, applies the frozen eligibility and
deduplication rules, reconstructs the inherited train/test source-row split,
and writes the six-class N/P cohort used by the current reports.

The source dataset is available from the [BCO-DMO Dataset 3228 landing
page](https://www.bco-dmo.org/dataset/3228). The page identifies the primary
file as `AllBottle.csv` and licenses the dataset under CC-BY-4.0. Please retain
that attribution when redistributing the processed subset. The raw file is
not committed by this project; download it from BCO-DMO and pass its path in
`params$raw_file` (or set `OCEAN_RAW_FILE`).

## 1. Required raw file and fields

Download **`AllBottle.csv` from BCO-DMO Dataset 3228** (“All bottle data …
ECOHAB-PNW project”). The preparation needs the following raw columns:

| Raw column(s) | Role in this cohort |
|:--|:--|
| `Cruise`, `Station`, `Survey`, `date`, `Event` | provenance retained in the output |
| `depth`, `depthID`, `bottle` | shallow-water eligibility, deduplication, provenance |
| `lon`, `lat` | retained location fields |
| `NO3_NO2` | observed response (y=\mathrm{NO_3+NO_2}) and inherited split |
| `H2PO4` | latent phosphate (u) and six ordinal classes (C) |
| `SiO2` | retained for provenance; **not** an eligibility requirement here |
| `temp1`, `temp2` | temperature, using `temp2` only as a valid fallback |
| `sal1_corrected`, `sal2_corrected`, `sal1_uncorrected`, `sal2_uncorrected` | salinity, preferring corrected sensors |

The code treats empty strings, `NA`, `NaN`, `nd`, and `ND` as missing. Numeric
conversion and the sensor-fallback rules below reproduce the existing Ocean
preparation scripts.

## 2. Parameters and paths

The default raw-file search includes the path used in the working project:
`real-data/raw/ocean/AllBottle.csv`.
For a clean checkout, place the downloaded file there or provide an explicit
path through `OCEAN_RAW_FILE` or `params$raw_file`. By default the output is
written to `real-data/ocean/`; set
`OCEAN_DATA_DIR` or `params$output_dir` to write elsewhere.


## 3. Read and normalize the raw BCO-DMO file

Temperature uses `temp1` first and valid `temp2` as fallback. Salinity uses
corrected sensors first, then uncorrected sensors. The valid ranges are the
same conservative screening ranges used by the original preparation code.


## 4. Eligibility, provenance, and deduplication

The primary cohort is restricted to shallow samples (`depth <= 10` m),
positive nitrate+nitrite, positive phosphate, and valid temperature/salinity.
Unlike an earlier sensitivity workflow, `SiO2` is not used as a filter.
Rows are deduplicated by `Cruise + depthID`; when `depthID` is missing, a
stable raw-row fallback key is used. The original `source_row` is retained so
the final split can be audited against the raw file.


## 5. Reconstruct the frozen 100/300 source-row split

The class-6 experiment uses a frozen 100/300 train/test cohort drawn from the positive-primary eligible pool. 
We select 20 training and 60 test observations within each prespecified nitrate stratum using 
`source_split_seed = 20260730`, thereby fixing the source rows independently of the final ordinal labels. 
We then define six phosphate classes on the full eligible pool and attach these labels to the selected observations. 
This separates cohort selection from the construction of the six-class ordinal input.


## 6. Construct the class-6 N/P role-swap cohort

For the ocean case, we define:

\[
y = \mathrm{NO_3+NO_2}, \qquad
u = \operatorname{standardize}\{\log(\mathrm{H_2PO_4})\}, \qquad
x = (\mathrm{temperature},\mathrm{salinity}).
\]

The six ordinal levels are constructed from five empirical phosphate-quintile cutpoints 
computed over the full positive-primary pool. The lowest quintile is further divided at 
standardized log-phosphate \(-1.5\). Both the cutpoints and the standardization parameters 
are estimated using the full pool, rather than only the 400 selected train/test observations.


## 7. Write the reusable cohort files

The two principal files are:

- `eligible_pool.csv`: all 1,393 positive-primary eligible records, with
  constructed six-class labels and both raw/log-standardized phosphate;
- `study1_data.csv`: the frozen 100/300 train/test cohort, calibration flags,
  and `sample_row_id`.

The metadata file records the raw-file checksum, source-row split seed,
calibration seed, cutpoints, standardization constants, and class counts.
This makes it possible to distinguish a changed raw download or preprocessing
rule from a changed model fit.


## 8. Reproducibility checks

The checks below are intentionally strict. In particular, class labels are
computed from phosphate, (y) is nitrate+nitrite, the calibration rows are
nested and selected without using (y) or exact (u), and the log transform
round-trips to the raw phosphate values.


The resulting `study1_data.csv` is the direct input expected by the Ocean
EIV-GP and baseline workflows. The model code should not recompute classes or
standardization from the 400 selected rows; it should read the generated
`meta.json` and use the stored `c`, `u_log_std`, `split`, and calibration flags.
