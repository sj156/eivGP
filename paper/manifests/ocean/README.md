# Ocean input manifests

- `inputs.csv`: row counts and SHA-256 hashes for the supplied 200/500 inputs,
  their calibration-coverage audits and the raw 1000-profile cohort.
- `upstream.json`: source provenance and workflow defaults.
- `validation.json`: scoped implementation checks and their run status.

The supplied 200/500 inputs use class-stratified nested calibration with
approximate Latin coverage of salinity and log pressure. Their local
`freeze_manifest.json` files record seeds and checksums. The input reader checks
this calibration contract before fitting. Preparation for the raw 1000 cohort
uses the same selector.

Paper sampling is four chains with 50,000 initial iterations per chain,
including 5,000 burn-in. Continuation preserves chain count and retained states.
The `cores` argument is a shared CPU budget. Default paper outputs are under
`results/<dataset>/`; `--output` selects a different run directory.

Input/provenance validation, a completed workflow and MCMC convergence are
separate checks. See the [Ocean guide](../../applications/README_ocean.md) for
commands and the [calibration design](../../applications/OCEAN_CALIBRATION_DESIGN.md)
for statistical scope.
