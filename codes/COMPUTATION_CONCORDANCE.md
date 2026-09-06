# Computation-to-code concordance

This file maps the posterior computation described in the manuscript to the
publication R implementation. Updated after the 2026-09-06 no-HMC revision;
the simulation manifest stores code hashes so that later runs can detect
changes.

## Common dense-GP target

- The finite GP vector is integrated out in every training transition.
- The covariance is
  `sigma2_eps * (I + rho^2 * R)` with either squared-exponential or Matérn
  ARD correlation on the joint `(X, U)` input.
- Likelihood evaluations use Cholesky solves and log determinants; the
  publication path adds no diagonal jitter and does not form an explicit
  inverse.
- `sigma2_eps` has the inverse-gamma full conditional with shape
  `2 + n / 2` and rate `0.05 + y' B^{-1} y / 2` on the standardized response
  scale.

Primary implementations:

- Study I: `00_study1_functions.R`, especially `gp_state_1d()`,
  `sample_sigma2_eps_1d()`, and `fit_eivgp_1d()`.
- Study II: `00_study2_functions.R`, especially
  `gp_loglik_integrated_general()`, `sample_sigma2_eps_general()`, and
  `fit_eivgp_ordprobit_fb()`.

## Study I sweep

The loop in `fit_eivgp_1d()` implements the manuscript order:

1. Draw `sigma2_eps` from its full conditional.
2. Draw the ordered deterministic thresholds from their feasible uniform
   conditionals.
3. Apply blocked Gaussian-prior elliptical-slice moves to missing `u`, with a
   periodic all-missing block.
4. Apply local slice moves in probability coordinates for well-resolved
   interior categories, or in direct u coordinates including the normal
   log prior for unbounded, extreme, or narrow categories. The coordinate
   choice depends on the category bounds, not the current u value.
5. Apply bounded componentwise slice moves to the GP log covariance
   parameters conditional on the current noise variance.
6. Check interval, calibration, and finite-state invariants throughout the
   sweep and before retention. Singleton index sets cannot expand to other rows.

This remains the default `noise_strategy="conditional"`. The optional
`noise_strategy="collapsed"` omits the initial noise draw, uses the analytic
inverse-gamma-integrated likelihood in every latent/kernel move, and regenerates
noise from the latest latent/kernel state before retention. The thresholds,
latent prior, kernel prior, calibration constraints, and posterior target do
not change. Dense and Schur calculations implement both strategies.

## Study II interwoven sweep

The `sampler_strategy = "interwoven"` loop in
`fit_eivgp_ordprobit_fb()` implements the manuscript order:

1. Draw truncated ordinal scores, cut points, and loading rows.
2. Integrate out `sigma2_eps` and apply score-conditioned local `U` moves and
   a parameter-only log-covariance move.
3. Periodically apply a score-conditioned joint move of all missing `U` and
   the GP log covariance parameters.
4. Periodically apply a score-marginal/prior-reference joint move of those
   quantities; redraw scores immediately.
5. Periodically update `(A, tau)` under the score-marginal ordinal likelihood;
   redraw scores immediately.
6. Assert that calibrated rows have never changed, draw `sigma2_eps` from its
   full conditional, check the state, and only then retain it. A transient
   violation cannot be repaired by overwriting calibrated rows at retention.

The regeneration order is mandatory for the partially collapsed chain. The
implementation fixes the ordinal-probit residual covariance to `Omega = I`;
it does not estimate residual proxy correlations. Triangular loadings together
with diagonal ARD are an explicit joint-model restriction. The public API now
requires an explicit loading choice with sparse/rank-deficient calibration.

Two further exact blocks are optional, with both frequencies zero by default:

- `joint_measurement_every`: each row's free loadings and all cutpoints are
  mapped from their original Gaussian/half-normal and ordered-uniform priors
  to independent standard normals. A small-row ESS uses the row's marginal
  ordinal likelihood. Scores are regenerated immediately afterwards.
- `loading_transport_every`: at fixed scores, missing latent rows are moved
  with proposed loadings while their standardized conditional latent residuals
  remain fixed. The Jacobian cancels the conditional Gaussian density. The
  residual comprises the integrated GP likelihood, missing-row marginal score
  densities, and calibrated-row conditional score densities. Calibration is
  untouched; noise is regenerated before retention.

These blocks are inserted after the scalar marginal measurement block and
before the noise draw. They require `sampler_strategy="interwoven"`. Control
validation rejects malformed schedules and frozen required blocks. No HMC,
gradient, new statistical prior, or likelihood approximation is introduced.

## Exact numerical savings

Both engines cache accepted GP determinant/quadratic-form summaries. Noise
draws do not require another factorization. At fixed kernel parameters, small
latent-block proposals can use the exact Schur complement of the unchanged
complementary covariance. Global/kernel moves remain dense.

The default selects this shortcut at n >= 120. Explicit TRUE/FALSE overrides
are `gp_block_schur` (Study I) and
`control_overrides=list(gp_use_block_schur=...)` (Study II). A NULL value
selects automatically. Study I additionally requires n >= 40 and small enough
blocks; Study II uses blocks smaller than n/2. Study I retries a failed
shortcut with the same dense covariance; Study II stops. Neither adds jitter.

Per-chain telemetry distinguishes full/complement/block calculations from
cache hits and residual calls. Effective samples per second—not proposal
acceptance or raw iteration counts—measure end-to-end efficiency.

## Posterior recovery and prediction

- `gp_integrated_mean_state_1d()` and
  `gp_integrated_mean_state_general()` implement the empirical-distribution
  GP linear-functional mean and covariance for `m(x,c)`.
- `sample_eiv_f_given_xu_fb()` and the corresponding Study I helpers recover
  `f(x,u)` by ordinary GP conditioning at retained states.
- Study II prospective `U | C` uses
  `sample_u_given_c_ordprobit_minimax()`: it samples the truncated marginal
  Gaussian score by minimax tilting and then samples the exact Gaussian
  conditional `U | S`. Solver, degeneracy, incomplete-sample, or unfamiliar
  dependency warnings stop the operation before it can be labeled exact.
  Finite Gibbs endpoints are diagnostic only.
- Predictive NLPD is computed from stored conditional Gaussian component means
  and variances by log-sum-exp, not by a kernel density estimate of predictive
  draws.
- Prospective prediction uses the fixed-fit convention; the interface does
  not update global posterior weights with a new unlabeled ordinal vector.
- Both engines fix the latent population law to N(0,I), independent of X.
  Public calibrated standardization is plug-in and is on by default; the
  low-level simulation engines preserve supplied latent units by default.

## Strict diagnostics

Both publication drivers now require four chains, Rhat <= 1.01, and bulk/tail
ESS >= 400 for required series. These thresholds are screening requirements,
not guarantees; MCSE and latent-integration sensitivity still need assessment.
Mean-specific MCSE and ESS now use `posterior::mcse_mean` and `ess_mean`;
bulk ESS is not substituted in the MCSE formula. Optional `mcse_limit` and
`mcse_sd_ratio_limit` arguments screen accuracy of the specified functional
mean. Their limits must be chosen for the application; they are not evidence
about model correctness or latent-integration bias.

Study I checks all free parameters and missing inputs, plus a separate panel
of m, predictive, and eligible f/U moments. Study II checks every ordinal
correlation and standardized cutpoint, the target panel in both calibration
regimes, and all free raw coordinates when calibrated. Conditional moments
use consecutive within-chain states and fixed common integration randomness,
not newly noised response draws or a short pooled reporting subsample.
Both target panels default to all retained reporting draws. Explicitly
truncated panels retain per-parameter incompleteness flags, including after
outer-list combination/subsetting, and cannot pass the full reporting gate.
Study II now monitors the integrated m conditional variance, and prospective
U moment diagnostics back-transform to supplied units when preprocessing was
used. Integration size and random-seed sensitivity remain separate checks.
Nonfinite diagnostics fail. Structurally fixed coordinates are excluded by
model definition, never because their computed diagnostics are inconvenient.

## Publication orchestration

- `simulation_helpers.R` contains the common frozen design and validation
  layer.
- `run_study1_simulation.R` and `run_study2_simulation.R` are the two master
  entry points.
- Publication parallelism is across replications; chains within a replication
  are serial to avoid nested forking.
- Public competitors are called through one adapter each in
  `03_study2_published_competitors.R`. Missing packages, inadmissible fits,
  diagnostic failures, exact-sampler underfill, and numerical failures are
  retained as failures rather than silently replaced.

## Verification status

The September 5 regression tests cover singleton/complete/no calibration,
tail support, narrow intervals, numerical failure handling, independent
conditional targets, exact dense/block agreement, API model contracts, and
strict diagnostic failures. Run `tests/test_study1_sampler_repairs.R`,
`tests/test_study2_sampler_repairs.R`, `tests/test_diagnostic_gates.R`,
`tests/test_study1_target_diagnostics.R`, and `tests/test_public_model_contract.R`.

`16_sampler_efficiency_benchmark.R` runs same-target paired timing checks and
saves the full fits. The fresh 3,000-iteration Study II checks still fail the
strict publication gates; publication readiness has not been established.

## Historical audit checks

The following focused checks passed on 2026-09-02:

- `09_experiment_design_validation.R`
- `13_study2_prospective_latent_validation.R`
- `14_mean_function_validation.R`
- `15_appendix_computation_validation.R`
- a reduced two-chain interwoven run of `05_study2_sampler_validation.R`,
  including its independent one-dimensional grid-target check

The full publication diagnostic gate still has to be satisfied by the final
workstation runs. Passing a smoke or reduced validation run is not evidence
for the numerical results reported in the manuscript.

The September 6 additions have independent quadrature, transformation/Jacobian,
calibration, regeneration-order, numerical-cache, and dense/Schur tests in
`tests/test_study1_noise_collapse.R` and
`tests/test_study2_joint_measurement.R`. Diagnostic precision is checked by
`tests/test_diagnostic_precision.R`. The literate package rebuild passes 178
tests and R CMD check with zero errors, warnings, or notes.

`17_no_hmc_mixing_pilot.R` runs both studies with frozen synthetic data, source
snapshots, four parallel chains by default, full-window target diagnostics,
and ESS/time comparisons. Its output is explicitly non-publication evidence.
The new pilots still fail strict convergence. Joint measurement rows appear
promising; the loading transport did not improve consistently and remains
experimental. See `../COMPUTATION-REVISION-20260906.md` for settings and results.
