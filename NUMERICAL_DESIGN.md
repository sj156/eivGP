# Numerical experiment design

Both studies use **100 training and 100 test observations per dataset**.
Publication uses **100 replications per setting**; development uses **3** to
exercise the workflow. Development defaults to every setting and the same
calibration grids, not a selectively smaller scientific design.

| Study | Setting | Calibration sizes |
|---|---|---|
| I | Balanced, eta=0 | 0, 20, 50 |
| I | Balanced, eta=1 | 0, 20, 50 |
| I | Imbalanced, eta=0 | 0, 20, 50 |
| I | Imbalanced, eta=1 | 0, 20, 50 |
| II | Primary q=2 | 50 |
| II | Primary q=4 | 0, 20, 50, 80 |
| II | Logistic misspecification q=4 | 50 |

Seven settings, 18 setting–calibration combinations: **1,800 publication
EIV-GP fits** or **54 development EIV-GP fits**, before other methods.
Each calibration set is nested within the same training dataset.
Study II has two latent dimensions; q counts ordinal proxies. Its primary
q=2 measurement mechanism is ordinal-probit, as assumed by the fitted model.
The imbalanced Study I generator retains its minimum of three training
observations in each ordinal category.

The existing default squared-exponential response kernel is unchanged.
The proposed pairwise-only kernel is not implemented or selected here.
Development retains four chains, 500 warmup plus 1,250 sampling iterations
per chain. **Publication MCMC budgets are not finalized**: existing numeric
defaults are placeholders, not an approved final publication specification.
No automatic extension or new prior/sampler change is made by this design update.

Designs live in codes/simulation_helpers.R; the synthetic-data R Markdown
documents and study masters use that configuration. The reusable eivGP 0.3.0
package deliberately excludes study design code and does not need rebuilding
for this change. Legacy standalone drivers are not the configuration entry point.

## Frozen-data transition

Existing frozen files, manifests, and figures are preserved. Older files may
have different sample sizes, calibration sets, or removed settings; they do
not become instances of the revised design through relabeling.
An incompatible existing manifest should stop data verification. Prepare
the revised frozen collection in a clean explicitly selected data directory
before running fits; do not delete or overwrite old files to force a match.
Development and publication still have separate default data locations:
shared storage/subsetting remains to be completed, although their generators,
seeds and per-replication designs now agree.
