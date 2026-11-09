# Study II disk-backed continuation

Use `tools/study2-continue-eivGP.R` for explicit continuation of saved fits.
`tools/study2-continue.R` is a compatibility alias. Both use the same
`tools/study2-disk-tools.R` implementation and require eivGP >= 0.4.0.9020.

These sensitivity experiments have separate output roots from the fixed-budget
paper replication. Original fits, frozen data and prior experiment libraries
remain read-only inputs. Generated artifacts are excluded from Git.

## Inputs and preparation

For an existing long pilot, supply `--pilot-root=/path/to/pilot`, or set
`EIVGP_STUDY2_PILOT_ROOT`. The portable default is
`results-sensitivity/study2-pilot-12cores`. A pilot root contains
`plan-pilot.csv` (including rep, calibration, fit and data columns) and committed
per-job `state.rds` files. Only saved states are selected; unstarted jobs are
not silently restarted.

Alternatively, use `--source=screen --screen-root=/path/to/sampler-screen` for
completed public ESS/HMC/HYBRID fits. The screen root contains `experiment.rds`
and each method's `fit.rds`; see [sampler-pilot instructions](study2-sampler-pilot-README.md).

Run commands from the repository root. Inspect the plan, then prepare a new
experiment with a pinned package library:

```sh
Rscript tools/study2-continue-eivGP.R --action=plan --pilot-root=/path/to/pilot
Rscript tools/study2-continue-eivGP.R --action=prepare --pilot-root=/path/to/pilot --library=/path/to/pinned/library
Rscript tools/study2-continue-eivGP.R --action=migrate
```

The default output is `results-sensitivity/study2-disk-continuation-v1`.
Override it with `--output=PATH` consistently in every command. For screen
inputs, also repeat `--source=screen` consistently and choose a separate output.
Preparation pins source/input hashes, R and dependency versions. Never replace a
library while a job uses it.

Migration is serial and reads one legacy fit at a time. The adapter attaches
missing public continuation metadata to supported internal Study II checkpoints;
data, priors, controls, retained draws, terminal states and RNG are preserved.
Sampling workers start only after import finishes.

## Sampling and restarting

After reviewing whether further sampling is warranted, request a bounded extension:

```sh
Rscript tools/study2-continue-eivGP.R --action=run --add=2000 --chunk=500 --round=1 --cores=12
Rscript tools/study2-continue-eivGP.R --action=report
```

Use `--jobs=job1,job2` to select prepared jobs. Repeating an identical request
resumes its fixed target without adding iterations twice. Another extension
requires the next round number after the previous round finishes.

Sampling reads the small committed index, terminal Markov/RNG states and bounded
new draw buffers. It does not load accumulated history or calculate R-hat/ESS.
HMC history is stored in disk blocks rather than the live terminal state.
Twelve cores allow three concurrent four-chain fits, with one BLAS thread each.
An atomic all-chain commit records progress; an interrupted uncommitted block is
recomputed from its saved RNG state. No thinning or repeated warmup is introduced.

## Assessment after sampling

Once sampling has stopped, assess histories separately:

```sh
Rscript tools/study2-continue-eivGP.R --action=diagnose
```

This serial phase reconstructs one retained history at a time and may use RAM
proportional to that history. It assesses scientific targets, ordinal measurement
invariants and applicable raw-coordinate checks using R-hat <= 1.01, bulk/tail
ESS >= 400 and MCSE/posterior-SD <= 0.05. Undefined required diagnostics fail.
There are no automatic ESS-based extensions in the sampling phase.

Persistent chain disagreement requires sampler/initialization review before more
iterations. Cheaper storage does not resolve the difficult rep 5, calibration 0
mixing issue. Retain diagnostic failures in comparisons.

Keep `fit.rds` and its adjacent `fit.rds.draws` directory together when copying
an archive. See `tests/application/test-study2-disk.R` for exact continuation
equivalence checks across ESS/HMC/HYBRID, calibration and score-retention settings.
