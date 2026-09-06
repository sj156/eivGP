# eivmixgp 0.2.0

- Added diagnose_eivgp(): full-window, target-aware R-hat, bulk/tail ESS,
  and estimator-specific mean MCSE with explicit failure classifications.
- Added continue_eivgp(): serializable terminal-state checkpoints, per-chain
  RNG restoration, frozen warm-up tuning, cumulative iteration schedules,
  and recomputed diagnostics after appending every new draw.
- Public fits now require no thinning (thin = 1). Multivariate defaults
  and preset retention settings also keep every post-warm-up draw.
- Diagnostic computation is shared by the package and simulation helpers.
- Added split/uninterrupted, serial/parallel, save/reload, repeated-resume,
  invalid-checkpoint and scientific-functional diagnostic regression tests.

# eivmixgp 0.1.0

- Added the stable `fit_eivgp()`, `predict_eivgp()`/`predict()`,
  `impute_eivgp()`, and `run_mixedgp_experiment()` interfaces.
- Added concise `print()` and calibration-aware `summary()` methods for
  fitted EIV-GP models.
- Made minimax-tilted accept--reject simulation, exact on successful
  completion, the default for prospective multivariate latent inputs; exact
  prior rejection and diagnostic finite Gibbs remain available explicitly.
- Required explicit multivariate latent dimension when calibration cannot
  supply it, and restricted raw-scale output to full-affine-rank calibration.
- Added squared-exponential and Matern kernels, multivariate numeric inputs,
  labeled ordinal-factor mappings, and exact-new-latent response prediction.
- Unified both numerical studies behind one runner with frozen-data collision
  checks, reproducible serial/fork RNG, audited competitor wrappers, and
  explicit sampler safety failures.
- Added generator/design fingerprints, deep truth validation, and atomic
  incremental manifests for frozen synthetic datasets.
- Added strict public count, logical, name, identification, and argument
  validation while preserving caller random-number streams.
- Kept plotting, table, and published-competitor packages optional so loading
  `eivmixgp` does not attach packages, modify the ggplot theme, or set thread
  environment variables.
