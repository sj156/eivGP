# Study II competitors: reproduce or resume

Run commands from the repository root. Install the locked dependencies and the
current package first:

```sh
Rscript tools/setup_numerical.R
```

After editing only this package, use `Rscript tools/setup_numerical.R --install-local`
instead of restoring dependencies again.

## Reproduce the complete comparison

The standard entrypoint includes the original UC-GP, LVGP and EzGP grid and the
calibrated q=4 comparisons at calibration sizes 20 and 50. Its default action is
a plan. Supply the same output/data roots on subsequent invocations:

```sh
Rscript tools/study2-competitors.R --action=plan --cores=12
Rscript tools/study2-competitors.R --action=run --cores=12
Rscript tools/study2-competitors.R --action=report --cores=12
```

Defaults are `~/eivGP-results/study2-competitors` for outputs and
`~/eivGP-results/data/study2` for frozen data. Override them with
`--output=/absolute/path` and `--data-root=/absolute/path`.
Rerun the same run command after an interruption. Completed scores are skipped;
saved predictions are scored directly; fitted models and committed measurement
checkpoints are reused. Measurement chains retain the published seeds and budgets.
Unfinished measurement work is checkpointed every ten minutes from startup by
default; checkpoint timing is operational and does not change random streams.
Prediction/scoring does not change the data, priors, fitted models, or integration
budgets. CRPS uses the exact mixture formula with bounded temporary arrays.
Scoring runs serially within a job while independent jobs share the core budget.

## Continue the interrupted local experiment without changing its originals

```sh
Rscript tools/study2-resume-competitors.R --action=plan --cores=12
Rscript tools/study2-resume-competitors.R --action=run --cores=12
```

This entrypoint finds `results/study2-competitors` in the checkout if present,
otherwise the standard output root above. It makes a verified copy at the source
path plus `-resumed`, then runs only the calibrated competitors by default. Use
`--methods=UC-GP,LVGP,EzGP,UC-GP-CC,LVGP-CC,EzGP-CC` to include the original grid.
The original experiment is never written. Repeating the resume command continues
its existing output; it does not copy the original experiment again.

Use `--source-root=/absolute/original` and `--output=/absolute/resumed` to select
other roots. Keep the same roots when restarting. `--data-root` selects the existing
frozen data collection. Plan starts no fits. Report operates on the resumed output.

Locks are recovered only after their local owners are verified inactive; lock
metadata is archived. Repeating the resume command while its verified local run
is active prints `Already running` and exits normally. Remote or unverifiable
owners still stop the command. The source-copy step checks
all recorded source locks and excludes their directories from the new copy.

The package explicitly permits the audited scheduling, scoring and cache-lifecycle
source changes when matching old protocols/checkpoints. It keeps the original
protocol and records the new runtime separately. Other source, dependency, data,
seed, iteration, prior, and integration changes remain subject to identity checks.
Real failures remain visible; warnings and errors are not suppressed.
