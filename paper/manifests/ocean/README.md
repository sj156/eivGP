# Ocean input and workflow identities

- `upstream.json`: source repository/commit and Ocean integration scope.
- `inputs.csv`: row counts and SHA256 for 33 input CSVs: 16 each for the 200/500 frozen designs and one raw 1000-profile cohort.
- `validation.json`: checks executed for this integration and their exact scope.

The package's installed `module-manifest.csv` records source module hashes. Output roots additionally store package/functions/data identities, per-task configs, checksums, completion receipts and fitting settings. The example HTML uses a fresh 200-profile single-fold quick run. The updated 500 and raw 1000 inputs have been checked; full-budget fits remain separate. The standalone fold script writes generated data under `reproduction/ocean/prepared/`, outside the versioned raw cohort. Original fitting archives are preserved outside this repository.
