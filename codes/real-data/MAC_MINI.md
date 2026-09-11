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
ADNI keeps four chains, using up to `min(N,4)` concurrent workers. With two workers,
all four chains are still run in batches. Folds and competitor fits run sequentially,
so budgets above four do not add parallel ADNI fits. During serial stages, some
of the budget will be idle. Start with two workers when sharing the Mac; use four
when it is available for the analysis. This is a worker budget, not an OS-wide
quota on unrelated applications.

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
writes `results/adni/repeat2/`; all three folds run sequentially. The explicit
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
