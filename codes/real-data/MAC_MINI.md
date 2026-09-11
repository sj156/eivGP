# Run the applications on a Mac mini

The Mac launcher uses the same notebooks and eivGP **0.3.1** as the Linux workflow.
It sets numerical-library threads before R starts and uses `caffeinate -i` while
R runs to prevent idle system sleep. It does not alter model, chain count,
warm-up, iterations, random seeds, or the frozen folds. It does not prevent
manual sleep, shutdown, or interruption of a running job.

## First check

In Terminal, enter your repository directory (adjust the path on your Mac mini):

```sh
cd "$HOME/Documents/GitHub/eivGP"
bash codes/real-data/run_mac.sh check --cores 4
```

The launcher detects ADNI data first in `real-data/adni/`, then directly in
`real-data/` (your current layout), then the older `codes/real-data/ADNI-toledo/data/`.
It reports missing packages/data and does not start MCMC. Ocean data can remain
unavailable while you run ADNI. To use another location, set `ADNI_DATA_DIR`.

If packages are missing or eivGP is not 0.3.1:

```sh
bash codes/real-data/run_mac.sh setup
bash codes/real-data/run_mac.sh check --cores 4
```

Setup needs internet and installs the tested eivGP commit plus the published
competitor packages. It may replace a different local 0.3.1 build. Finish setup
before starting a run, and keep that package installation stable when resuming.
Setup session information is saved in `results/setup/`.

If R is not installed, use the [official R macOS installer](https://cran.r-project.org/bin/macosx/)
matching your Mac's architecture and macOS version (Apple silicon or Intel).
Use `RSCRIPT_BIN=/full/path/to/Rscript` if it is not on PATH. RStudio and Pandoc
are unnecessary for this launcher. If installation reports a source compilation
failure, follow the compiler guidance for your R release on the official R page;
installing a matching binary package avoids compiling it locally.

## One setting for the CPU budget

```sh
bash codes/real-data/run_mac.sh adni --cores 4
```

`--cores N` is the total analysis-worker budget. The default is 2; the option
also accepts `--cores=4`. There is no need to set separate chain/thread variables.
The launcher limits OpenMP, OpenBLAS, MKL, and Accelerate numerical threads to one.
ADNI keeps four chains per fit. The scheduler divides the budget across the
selected repeat/fold jobs, maximizing concurrent chain workers and preferring
fewer simultaneous folds on ties to reduce memory use. For one repeat, 12 cores
run three folds with four workers each. Select both repeats with `--repeats 2,3`
to expose six jobs: 16 cores run four folds at a time; 56 cores run all six folds
with up to 24 chain workers. The remainder is idle because no additional fits
are requested. During competitor/reporting stages, usage may be lower. This is
an analysis-worker budget, not an OS-wide quota.

Each concurrent fold holds a separate fit in memory. If memory pressure or swap
becomes substantial, limit concurrent folds with `ADNI_MAX_FOLD_WORKERS=2` (or 1),
while keeping `--cores` as the total budget. No automatic hardware-memory estimate
is assumed. Do not increase chains or repeats merely to occupy spare cores.

The R runner accepts the same option on either platform:

```sh
Rscript codes/real-data/run_application.R adni --cores 4
```

Use the shell launcher on macOS to set numerical-library limits before R starts
and keep the Mac awake. Advanced automation can set `EIVGP_CORES`; `--cores`
overrides it. Existing direct-R commands using `EIVGP_N_CORES` remain supported
when neither the unified variable nor the flag is supplied.

## Smoke test, then production

```sh
bash codes/real-data/run_mac.sh smoke --cores 2
bash codes/real-data/run_mac.sh adni --cores 4
```

Smoke writes `results/adni/smoke_repeat2/`. Production defaults to repeat 2 and
writes `results/adni/repeat2/`; all three folds can run concurrently within the total core budget. The explicit
`adni` action selects production even if EIVGP_SMOKE_TEST was left in your shell.
The smoke test uses tiny MCMC and optimizer budgets: numerical competitor
failures must be inspected but are not evidence of production performance.

For a run that continues after Terminal closes:

```sh
mkdir -p results/logs
nohup bash codes/real-data/run_mac.sh adni --cores 4 > results/logs/adni-repeat2.log 2>&1 &
echo $! > results/logs/adni-repeat2.pid
```

Watch progress in another Terminal:

```sh
tail -f results/logs/adni-repeat2.log
```

Ctrl-C exits `tail` only. For a foreground analysis it interrupts the analysis.
For a background job, use Activity Monitor to identify and stop the analysis's
R process and its workers. The repeat's `.run-lock/owner` records the main R PID.
Confirm those processes have stopped before restarting. An interrupted segment after the
latest checkpoint may be lost; normal checkpoint resumption retains completed
sampling segments. Do not start two processes on the same repeat output folder.

After repeat 2 finishes:

```sh
ADNI_REPEAT_ID=3 bash codes/real-data/run_mac.sh adni --cores 4
bash codes/real-data/run_mac.sh report
```

The final command combines completed repeats 2 and 3 without fitting models.
See [ADNI-toledo/README.md](ADNI-toledo/README.md) for report contents and method
failure handling. Main figures/tables are in each repeat's `main/`; detailed
results are in `appendix/`. Combined tables are in `results/adni/combined/`.

## Resume or move a run

Rerun the same command to resume compatible checkpoints. The core budget can
change without changing the statistical configuration. To resume work from
another machine, transfer the whole repeat directory, including run_config.csv,
fit_latest.rds, and case-study caches. Keep the installed package code compatible;
package fingerprints can differ across installations and reject comparator-cache
reuse. Do not copy a live, changing checkpoint directory.

If a killed job leaves `.run-lock`, confirm the old R process and its workers
have stopped before removing that repeat's lock directory. A completed workflow
is not necessarily a converged or complete four-method comparison: inspect
method status and diagnostics.

## Ocean

```sh
OCEAN_DATA_DIR=/path/to/prepared-ocean bash codes/real-data/run_mac.sh ocean --cores 4
OCEAN_DATA_DIR=/path/to/prepared-ocean OCEAN_QUICK=false bash codes/real-data/run_mac.sh ocean --cores 4
```

The default ocean mode is quick. It still needs `study1_data.csv` and `meta.json`.
Outputs default to `results/ocean/`. `ADNI_OUTPUT_DIR` and `OCEAN_OUTPUT_DIR` can
point to another local drive. Private data under `real-data/` and generated
outputs under `results/` are excluded from Git; Git does not transfer them.

## Parallel fold scheduling

Use one coordinator for the requested repeats:

```sh
# Inspect allocation without reading data or fitting:
Rscript codes/real-data/run_application.R adni --cores 16 --repeats 2,3 --plan
# Fit both repeats in a shared job pool:
Rscript codes/real-data/run_application.R adni --cores 16 --repeats 2,3
# A larger machine uses the same interface:
Rscript codes/real-data/run_application.R adni --cores 56 --repeats 2,3
```

`--repeats` defaults to ADNI_REPEATS, then ADNI_REPEAT_ID, then repeat 2. A
single repeat has three fold jobs (at most 12 chain workers); both repeats have
six (at most 24). For both repeats: 12 cores -> 3 folds x 4 workers; 16 -> 4 x 4;
56 -> 6 x 4. Folds are dynamically dispatched as slots become free. Statistical
chain counts, seeds, sampler controls, and iteration budgets are unchanged.
The Mac shell launcher also accepts `--repeats 2,3`.

The coordinator holds each selected repeat's `.run-lock`. Each worker holds its
own `.fold-lock` and writes only that fold's files. Follow `fold_N/worker.log`,
`PROGRESS.md`, or `CURRENT_STATUS.txt`; the terminal need not print every worker
message. The coordinator writes combined reports only after requested workers
finish. `scheduler-repeat2-3-status.csv` in the output base records job status.
A failure preserves other completed checkpoints, stops combined reporting, and
records its error; competitor failures handled within a fold remain explicitly
reported by the existing case-study policy.

An already-running process retains the code it loaded. Do not start the new
coordinator alongside it. To switch, allow a checkpoint to complete, interrupt
the old run, confirm its R processes/workers have stopped, and rerun with the
new command. Remove a stale repeat/fold lock only after confirming its owner is
gone. Existing compatible checkpoints are reused; an interrupted unsaved segment
must be recomputed. Changing concurrency does not change checkpoint compatibility.

Direct Rmd rendering supports automatic fold scheduling within its selected
single repeat via EIVGP_CORES; the R runner coordinates multiple repeats.
For an all-fold smoke check only, set EIVGP_SMOKE_TEST=1 and
ADNI_SMOKE_ALL_FOLDS=1; smoke outputs remain separate from production.

## Fresh start with all three repeats

Repeats 1, 2, and 3 are supported. All three give nine fold fits per method.
A 12-core budget runs three folds concurrently with four chain workers each;
a 56-core budget can run all nine folds with up to 36 chain workers. This does
not create additional independent participants or remove earlier model-selection
history associated with any of the frozen repeats.

Stop any old run and confirm its workers have exited first. From the repository
root, the following deletes prior ADNI outputs from both standard locations,
including checkpoints and comparator caches, then starts all three repeats:

```sh
export ADNI_OUTPUT_DIR="$PWD/results/adni"
rm -rf -- "$ADNI_OUTPUT_DIR" "$PWD/codes/real-data/ADNI-toledo/outputs"
Rscript codes/real-data/run_application.R adni --cores 12 --repeats 1,2,3 --plan
EIVGP_SMOKE_TEST=0 Rscript codes/real-data/run_application.R adni --cores 12 --repeats 1,2,3
# After all fits finish, using the same ADNI_OUTPUT_DIR:
Rscript codes/real-data/ADNI-toledo/report_adni.R --repeats 1,2,3
```

Delete outputs only for an intentional fresh restart. To resume an interrupted
run, keep the same ADNI_OUTPUT_DIR and rerun without the deletion command.
Custom output directories from earlier runs must be removed separately.
The `--repeats` flag selects CV assignments; by itself it does not reset existing
runs. The default repeat remains 2 for backward compatibility.
