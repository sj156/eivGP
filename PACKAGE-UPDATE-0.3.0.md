# eivGP 0.3.0

- Restores the package name `eivGP` and the single installable repository
  directory `eivGP/`. Canonical code remains in `codes/`; the active build is
  `Rscript --vanilla litr/render-package.R`.
- Replaces the independent signal/noise prior with common variance
  `V ~ IG(3,2)` and signal fraction `r ~ Beta(32,8)`. Integrates out `V` while
  sampling, then recovers it conditionally for joint posterior prediction.
- Uses continuous numeric-input kernel coefficients and a finite dictionary
  for latent-input kernel coefficients. Updates dictionary coordinates from
  their categorical conditional distributions.
- Places Dirichlet priors on ordinal category probabilities, inducing the
  threshold prior. Elliptical slice and scheduled transport updates target
  this posterior without HMC or thinning.
- Retains configurable fixed iteration budgets, target-aware diagnostic
  warnings, and explicit compatible-checkpoint continuation. General defaults
  are four chains, 500 warmup and 1,250 retained transitions per chain.
- Keeps numerical-study generators, runners, budgets, frozen data, and
  reporting scripts in the repository experiment layer, outside the installed
  package API. Optional published-competitor adapters remain reusable.

This is a posterior change. Earlier `eivGP` and `eivmixgp` fits, caches, and
checkpoints must be refitted for 0.3.0; changing their labels is not sufficient.
Existing study results are not recomputed by this release. Short verification
runs are software checks, not evidence that the scientific diagnostics pass.
