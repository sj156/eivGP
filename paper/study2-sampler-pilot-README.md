# Study II sampler screen: replication 5, calibration 0

`tools/study2-sampler-pilot.R` compares fresh ESS, HMC and HYBRID fits on the same
archived data, priors and four initial chain states. Defaults are **4,000 total
iterations and 2,000 warmup iterations** per chain. This is a sensitivity screen;
the fixed-budget paper results are kept separately.

## Prepare reference inputs

The screen needs a `plan-pilot.csv` with rep, calibration, fit and data columns.
The follow-up planner can create it from an existing Study II run:

```sh
Rscript tools/study2-followup-plan.R --source=/path/to/saved/core-run --data-root=/path/to/frozen/study2 --output=results-sensitivity/study2-pilot-12cores
```

The planner starts no fits. Its pilot selection covers replications 5, 8, 16, 25
and 42 at calibration 0 and 20; the sampler screen selects only rep 5, calibration 0.
Reference fits and frozen datasets must already exist and match.

Use `--existing-pilot=PATH` or `EIVGP_STUDY2_PILOT_ROOT` for another manifest root.
The portable default is `results-sensitivity/study2-pilot-12cores`.

## Run the matched screen

From the repository root:

```sh
Rscript tools/study2-sampler-pilot.R --action=plan
Rscript tools/study2-sampler-pilot.R --action=prepare
Rscript tools/study2-sampler-pilot.R --action=validate
Rscript tools/study2-sampler-pilot.R --action=run --cores=12
Rscript tools/study2-sampler-pilot.R --action=report
```

The default output is `results-sensitivity/study2-sampler-screen-4000-v1`.
Preparation requires a fresh output root, snapshots the package and installs an
isolated library. Reuse the same output and budget options for validation, run
and reporting. Do not prepare again in an existing experiment or replace a
library used by an active job.

Runtime validation checks the pinned source, package/dependency identities,
archived input hashes, normalized priors and initial states. All retained draws
are saved. HMC adapts only during warmup; HYBRID adds a full ESS sweep every
20 iterations. Other archived model/sampler settings are preserved.

Twelve cores allow three concurrent four-chain fits. Per-chain checkpoints are
saved every 2,000 iterations or 900 seconds. Repeating the same run command resumes
incomplete checkpoints and reuses completed fits. After verifying processes have
stopped, `--recover-stale-locks=true` allows recovery of known stale locks.

A longer screen requires a fresh output root and explicit `--iterations` and
`--warmup` values, repeated consistently across preparation and execution.

## Diagnostics and outputs

Each method writes its fit, scientific/invariant/raw diagnostics, target and
A[1,1] traces, timing and applicable HMC diagnostics. Reporting writes
`sampler-comparison.csv` and, when assessments are available,
`scientific-target-agreement.csv`.

Required scientific/invariant diagnostics use R-hat <= 1.01, bulk/tail ESS >= 400
and MCSE/posterior-SD <= 0.05. Raw coordinates are reported but are not a
calibration-zero acceptance gate. HMC/HYBRID additionally require no flagged
divergences and valid finite step sizes/acceptance probabilities. Inspect traces
before interpreting a pass.

Minimum bulk ESS per fitting second is ranked only among passing fits with
comparable uninterrupted timing; warmup cost is included. Differences in posterior
means relative to combined MCSE remain descriptive when either method fails.
Missing assessments are never counted as passes.

## Recorded pilot outcome

The local 4,000/2,000 screen completed for all three samplers. All failed the
diagnostic gates, and the boundary configuration persisted in chain 4. This
single dataset and seed block therefore supplies no passing sampler to select.
The earlier 80,000/160,000-iteration continuation states remain separate
historical evidence, not matched timing controls.

Generated results, installed libraries and source snapshots are excluded from
Git. A fresh checkout must create its own inputs and prepared experiments.
For disk-backed extensions after diagnostic review, see
[disk-continuation instructions](study2-disk-continuation-README.md).
