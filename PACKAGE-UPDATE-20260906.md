# eivmixgp 0.2.0: diagnostic and continuation workflow

The EIV-GP implementation is distributed under the R package name `eivmixgp`.
The new archive is `package-release/eivmixgp_0.2.0.tar.gz`; it supersedes the
older 0.1.0 archive, which is preserved rather than overwritten. This is a
local source release, not a CRAN or GitHub publication.

## Changes

- `fit_eivgp()` retains every post-warm-up draw; `thin > 1` is rejected.
  Multivariate engine defaults and preset retention settings are also one.
  Low-level historical thinning support remains for reproducing older analyses;
  the current public workflow and publication master defaults do not thin.
- `diagnose_eivgp()` checks the complete retained window. It reports free
  parameters, missing training latent inputs and their squares, measurement
  invariants, and scientific target signatures on a user-selected panel.
  Scientific signatures include conditional mean/GP variance for m(x,c),
  predictive second moments, prospective latent moments when calibrated, and
  f(x,u) conditional mean/variance when an eligible calibrated U panel is supplied.
- User functionals such as scientific contrasts can be supplied as chain-aligned
  `additional_series`. The report uses rank-normalized split R-hat, bulk/tail
  ESS, and mean-specific MCSE, with explicit per-functional failure reasons.
- `continue_eivgp()` resumes the terminal state of each existing chain, including
  ordinal scores even when score histories were not requested. It preserves RNG
  kinds/streams, frozen warm-up tuning, accepted likelihood caches and absolute
  update schedules; it appends draws within chains and recomputes raw diagnostics.
- Checkpoints survive RDS serialization and reject modified data/draws/settings.
  No seed reset or changed model options are accepted by the continuation API.
- Canonical diagnostics now live in `codes/00_diagnostics.R`, shared by the
  package and simulation helpers. Public workflow code is in
  `codes/00_mcmc_workflow.R`. The literate build copies these exact modules.

## Intended use

```r
library(eivmixgp)
# fit must have been created with the new fit_eivgp(), preferably >= 4 chains.
report <- diagnose_eivgp(fit, X = X_panel, C = C_panel)
report
report$table
saveRDS(fit, "fit-checkpoint.rds")
# If investigation supports further sampling:
fit <- continue_eivgp(readRDS("fit-checkpoint.rds"), n_iter = 2000L)
report <- diagnose_eivgp(fit, X = X_panel, C = C_panel)
```

`n_iter` here means additional transitions per chain, all retained. There is
no second warm-up. Fits from earlier releases lack terminal checkpoints and
must be refitted; discarded historical draws cannot be recovered. Use the same
package version and numerical environment for reproducible continuation.
Checkpoints are returned when a fit finishes; periodic crash recovery inside
a running fit is not implemented.

## Interpretation and limits

The default screen requires at least four chains, R-hat <= 1.01, bulk and tail
ESS >= 400, and mean MCSE/posterior SD <= 0.05 for monitored nonconstant
functionals. These thresholds are screening defaults, not high-stakes accuracy
guarantees. Raw-coordinate and invariant results are kept distinguishable.

- Poor exploration: inspect traces, chain-specific targets, calibration and
  symmetries; running longer alone may not help.
- Insufficient precision: chain agreement passes the R-hat screen but ESS or
  mean MCSE is inadequate; extension followed by rechecking is reasonable.
- Insufficient diagnostics: no pass can be issued for too few chains, invalid
  values, or unresolved constant functionals.
- Screen passed: only the monitored functionals meet the chosen thresholds.
  This does not establish convergence, identification or model adequacy.

Common integration randomness avoids adding fresh prediction noise to the
mixing screen. Its finite numerical-integration error is not part of the MCMC
MCSE: check larger `n_latent` and different integration seeds. Quantile MCSE,
simultaneous bands and decision-specific absolute accuracy require additional
assessment. The new workflow does not establish good mixing for all datasets.
No theory or manuscript sections were changed in this package update.

## Verification

On macOS arm64, R 4.6.0:

- Final source-package `R CMD check --no-manual`: 0 errors, 0 warnings, 0 notes.
- Packaged tests: 257 passed assertions, 0 failures/warnings/skips.
- Exact split/uninterrupted equality for both Study I noise strategies and
  Study II legacy/interwoven strategies, including optional joint measurement
  and loading-transport schedules, SE/Matérn kernels and multiple X/U inputs.
- Exact equality after save/reload, serial/parallel continuation, repeated
  extensions, and changes in the caller's RNG kind; caller RNG is restored.
- Warm-up adaptation frozen, no thinning, unchanged calibration, checkpoint
  corruption rejection, and complete-window target/MCSE checks.
- Existing diagnostic gates, diagnostic precision, Study I target diagnostics,
  Study I noise-collapse/quadrature and Study II joint-update/Jacobian tests pass.

Windows/Linux installation and cross-platform numerical reproducibility have
not been tested in this update. These are implementation tests, not new
scientific simulation results or evidence that a particular fit has converged.

Archive SHA-256:
`3e82932464bbe34998de9a6c33edb23e90278d0b19045f7fe77ab98725bc8eec`.
