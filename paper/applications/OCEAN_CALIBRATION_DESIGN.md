# Ocean calibration design

## Scope and rule

The supplied 200/500 five-fold inputs use the selector described below. The
cohort, train/test partitions, measured records and test-U availability are
fixed. Preparation computes nested training-calibration indicators from recorded
class and X coordinates. It does not select a new cohort or screen model fits.

Within each recorded C class, allocate 20%/50% of training profiles. For the
current 500 cohort this means 20/50 of 100 training profiles per class, or 80/200
per fold. Use only salinity and log pressure to choose profiles. At each level k,
map each X axis to k empirical equal-probability rank bins within that class;
ties receive their average rank. Minimize the sum of squared differences between
bin counts and one over both axes and both calibration levels. Membership labels
(hidden / O50 only / O20+O50) enforce exact quotas and nesting throughout.

Three fixed-budget seeded searches with membership exchanges and a final
improving-exchange pass optimize X coverage only. Each search uses
max(2000, 30 * class size) proposals. Early exit occurs only at zero occupancy
loss. Seeds are fixed before fitting and are not selected using U, Y, convergence
or predictive/imputation scores. Stable row identifiers determine reproducible
ordering, not a physical-distance criterion. This heuristic does not prove a
global optimum. Selection uses recorded C/X without oxygen-derived groups.

## Why this is not guaranteed strict LHS

Classical LHS generates coordinates with one point in every marginal probability
stratum. Here points must be existing profiles; their X coordinates are correlated
and sometimes tied, and O20 must remain a subset of O50. A feasible joint nested
Latin subset is not guaranteed, and the bounded heuristic may also miss a feasible
exact solution. We keep equal recorded values together rather than manufacturing
coverage by splitting ties. Empirical-quantile coverage is relative to the observed
class-conditional X distribution; it is not equal-width spacing or a uniform
joint distribution over a rectangular domain. Marginal coverage also does not
ensure all joint X regions or rare subgroups are covered.

The accurate label is **C-stratified nested empirical-rank X coverage** or
**approximate Latin coverage**. `calibration_X_coverage.csv` reports actual bin
occupancy and marks a margin exact only if every bin has one selected profile.
No strict-LHS claim is made when that check fails. The method's extra complexity
should be justified by a coverage objective; proportional class-random sampling
remains a simpler valid alternative.

## Alignment with the manuscript

The 26 October draft's posterior-update condition F2 allows a parameter-free
calibration mask based on recorded X and C, including dependent fixed-quota
masks. Its selection discussion around equation (17) distinguishes this from
using unrevealed within-category U. Conditional on the supplied cohort/design,
the new mask is a function of X, C, fixed identifiers and independent random
numbers, and does not depend on missing U or model parameters. Consequently the
mask probability can be omitted from the missing-input likelihood under those
conditions. Optimization of observed-X coverage does not change that argument.
The mask is selected without unrevealed within-category input values.

The draft's synthetic LHS generates X; it is not its calibration algorithm.
Study I uses independent random calibration prefixes; Study II uses randomized
(C1,C2) strata. The Ocean rule follows the same observable-information
restriction but is not an exact replication of either simulation design.

The draft's consistency theorem uses a separate iid experiment and calibration
positivity condition. It does not automatically apply to this finite, balanced,
fixed-quota design. Better X coverage does not validate the Gaussian input law,
the measurement model, individual U posterior calibration, or physical-surface
identification over unobserved regions.

## Remaining design concerns

- The supplied cohort and folds were selected/balanced using exact oxygen. This
  design addresses the calibration-mask mechanism conditional on those records;
  it does not establish a representative ocean population or external transfer.
- Source metadata records development screening. Retain that provenance;
  do not describe the frozen cohort as selected without prior exploration.
- Changing calibration requires refitting EIV-GP and CC comparators. The original
  categorical fits use unchanged X/C/Y, but result reuse must pass the workflow's
  provenance checks; do not overwrite old identities to force reuse.
- Default commands use the supplied calibration inputs. Custom prepared inputs
  use the same explicit input/output paths for every method and report stage.
- The preparation and coverage checks validate an implementation, not statistical
  performance or MCMC convergence. Mask sensitivity can be considered separately;
  a seed must not be chosen because it produces better model scores.
