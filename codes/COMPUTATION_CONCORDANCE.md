# eivGP 0.3.1: computation-to-code concordance

## Transition revision (posterior unchanged)

For deterministic thresholds, let P_k = Phi(tau_k), P_0=0, P_m=1.
At fixed U and neighboring cutoffs, the split
T_k=(P_k-P_{k-1})/(P_{k+1}-P_{k-1}) has a Beta(alpha_k,alpha_{k+1})
prior conditional. Its full conditional is this distribution restricted by
max(U_i : C_i <= k) <= tau_k < min(U_i : C_i > k), intersected with the
neighbor interval. Empty groups impose no additional bound. Gaussian prior
factors for U and the collapsed GP likelihood are constant in this update.
One sequential sweep updates every internal cutoff, using inverse Beta CDFs
and both log tails, then maps probabilities back to Gaussian stick coordinates.
This is a coordinate Gibbs kernel for the existing posterior, not a new prior.
The ordinary latent block and joint cutoff/input transport are retained.
`threshold_update="ess"` selects the old cutoff step for comparisons.

For the probit joint move, a subject block and one measurement row change.
All other U and all GP kernel parameters stay fixed. The same exact Schur
calculation used by ordinary U blocks is therefore valid. The ordinal residual
is item j for all subjects plus other items for the selected subjects; omitted
terms are constant along the slice. The old dense fallback remains active.

New fits/checkpoints identify sampler 0.3.1; 0.3.0 checkpoints are not silently
reinterpreted. The prior/target specification below remains unchanged.

The reusable sampler lives in 00_sampler_v030.R; 00_sampler_v030_api.R
validates/scales data and packages posterior draws. The two historical
fit_eivgp_1d/fit_eivgp_ordprobit_fb entry points forward to this implementation.
The installed package is eivGP/. Numerical designs, competitors selected for
a study, data generation, replications, and reporting remain repository scripts.

## Target and transitions

- Signal covariance is V*r*R and response-noise variance is V*(1-r).
  Independently, V ~ IG(3,2) and r ~ Beta(32,8), on the standardized response
  scale. IG(a,b) has density proportional to v^(-a-1)*exp(-b/v).
  r is a pointwise variance fraction, not realized regression R-squared.
- GP values, V, and ordinal-probit scores are marginalized in the primary
  chain. With B=(1-r)I+rR and Q=y' solve(B,y), the log likelihood, up to a
  parameter-independent constant, is
  -0.5*logdet(B) - (a+n/2)*log(b+Q/2).
- The X-kernel log inverse squared-distance coefficients have independent
  N(log(.5),1.5^2) priors. Only U-associated coefficients have a finite
  dictionary. U itself stays continuous. SE and Matérn use the same ARD
  distance convention across fitting and prediction.
- The default U dictionary has five equal-weight quantiles per coordinate;
  coordinate categorical Gibbs avoids enumerating the Cartesian product.
  A user-specified matrix defines complete parameter vectors instead.
- Category probabilities have Dirichlet(2,...,2) reference priors. Gaussian
  stick coordinates generate these priors. Deterministic thresholds are
  normal quantiles; probit thresholds multiply them by sqrt(1+||A_j||^2).
  Free loadings have N(0,2.5^2) priors, with half-normal positive diagonals
  when the declared structural rule requires them.
- Every sweep visits all missing U coordinates in state-independent random
  subject blocks (at most eight subjects by default). Threshold models use
  interval-preserving Gaussian coordinates with the GP likelihood as ESS
  residual. Probit models use the population Gaussian reference and the
  GP likelihood times the score-marginal ordinal likelihood.
- Probit measurement parameters update jointly within each item; threshold
  cutoffs use the coordinate Gibbs step above. At fixed U neither update
  requires a GP factorization. Every ten sweeps, the reference schedule adds
  coordinated measurement/input moves. Threshold transport moves every
  missing input, includes category factors only for missing subjects, and
  preserves calibrated inputs. Probit cross-block moves combine one
  cyclically selected measurement row with a random subject block.
- Continuous X coefficients and the signal fraction update jointly by ESS,
  then the dictionary index updates categorically. There is no random-walk
  acceptance-rate adaptation and no HMC.
- Optional dictionary_mode="marginal" uses a sum of likelihoods throughout
  the continuous blocks and redraws the full index before recovery. It never
  averages covariance matrices. max_dictionary_size guards exponential cost.
- Every post-warmup state is retained. V is recovered from
  IG(a+n/2,b+Q/2), followed deterministically by sigma2=V*(1-r).
  Optional scores are recovered only for storage, never used by the next sweep.

All priors and sampler controls are explicit fit_eivgp arguments. Reference
choices need scientific prior-predictive and sensitivity assessment; exact
invariance is not evidence that mixing will be good in every application.

## Exact calculations and posterior outputs

Dense evaluation uses Cholesky solves without an explicit inverse or jitter.
Within a latent block the unchanged complementary factor supports exact Schur
calculations. Setup costs O((n-b)^3); candidates cost O(n^2*b+n*b^2+b^3).
Invalid shortcuts retry the unchanged dense covariance. Stale complementary
factors are rejected using their input and kernel keys. Global parameter
moves and threshold transports generally require dense evaluation.

Existing predictive helpers are reused by the exact correspondence
rho^2=r/(1-r), sigma2=V*(1-r). The retained native V,r,J,pi draws accompany
these compatibility coordinates. Empirical integration of m(x,c) averages
both GP means and covariance, rather than means alone. Prospective U uses
the existing exact truncated-Gaussian construction; fixed-fit prediction
does not reweight the training posterior using new ordinal inputs.

No probability clipping or adaptive jitter is used. An out-of-support
candidate is an ordinary ESS rejection. Undefined calculations, failed exact
draws, and exhausted ESS searches are explicit errors, not convergence gates.

## Diagnostics and continuation

Report rank-normalized R-hat, bulk/tail ESS, mean MCSE, and target-specific
diagnostics; native variance, signal fraction, category probabilities, and
dictionary occupancy remain visible. Unvisited nontrivial dictionary states
are not silently excluded. Movement angles, likelihood evaluations, and
full/block factorization counts accompany each chain.

Independent chains use reproducible serial/fork execution. Checkpoint schema
2 records the 0.3.1 sampler version, terminal collapsed state, RNG state, and
valid covariance cache. Public continuation checks the full fit signature,
freezes model/control choices, and appends within-chain draws. Pre-0.3.0
posteriors/checkpoints must be refitted, not relabeled.

Diagnostic warnings do not prevent completion of a requested finite budget.
All stored draws are retained; postprocessing integration budgets are separate.

## Verification and next-stage experiments

codes/tests/test_v030_core.R independently checks likelihood/scale integration,
prior transforms, dictionaries, transport factors, calibration support,
prediction covariance, and block/dense equivalence including failed shortcuts.
codes/tests/test_v030_workflow.R checks the public interface, all-draw retention,
stored-score recovery, checkpoint tampering, and exact resumed/uninterrupted
and serial/forked reproducibility.

Short verification fits are software tests, not numerical-study results or
proof of adequate scientific mixing. Experiment drivers use the new fitting
entry points but their study designs and result values are unchanged. Existing
cached fits and historical ablations are not evidence for the new posterior;
the experiment revision should audit their provenance, prior matching, and
reporting fields before starting new study runs.
