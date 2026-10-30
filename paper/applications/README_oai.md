# OAI reproduction and storage guide

The OAI application has six command-line and helper files for data preparation,
EIV-GP fitting, competitor fitting and report generation. `tools/oai-tools.R`
defines the application functions for inputs, execution and reporting. The
eivGP package supplies the samplers, numerical prediction functions, native
comparator implementations and shared CPU-pool utilities.

## 1. Install and supply inputs

Follow the repository installation instructions, reinstalling this checkout.
Place the privately supplied 46 files directly in `data/oai/`; see
[the data contract](../../data/oai/README.md). An external directory works with
`--data-root=/local/oai`. Paths in the code are relative to the checkout; commands
can be invoked from another working directory.

The default scientific design uses baseline quantitative MRI cartilage thickness
as U, fixed-location medial radiographic joint-space width as Y, continuous
age/BMI/height as X, and eight native central-reading radiographic grades as C.
Only the documented adjacent lateral-femoral-sclerosis grades are pooled in the
supplied cohort. Metadata provides all coordinate labels, grade maps and
observation design; another supplied p/q/d configuration does not require
rewriting the fitting code.

## 2. Paper execution

```sh
Rscript tools/oai-eivGP.R --plan=true --profile=paper --cores=8
Rscript tools/oai-eivGP.R --profile=paper --cores=8
Rscript tools/oai-competitors.R --profile=paper --cores=8
Rscript tools/oai-report.R --profile=paper --cores=8
```

The EIV-GP and competitor commands may run concurrently in two Terminals against
the same output root. Run the final report command after both finish. `--cores`
is the sole resource control, is per command, and is capped at detected capacity. EIV-GP fits use
min(cores,4) chains; above four cores, multiples of four permit complete resident
EIV groups. Each native competitor job uses one core. Lower the budget if memory
is limited: two commands each requesting eight cores need a combined budget of
16 workers. Changing only resources for saved-only reporting does not refit.

The two fitting groups own separate locks, task folders, logs, method reports
and progress summaries. A short initialization lock protects the common frozen
inputs, manifest and fitting identity. Duplicate commands for the same group are
blocked. The combined report temporarily requires both groups to be idle.

Console messages identify `START`, `DONE`, `REUSED` and `FAILED` jobs with fold,
calibration, method and elapsed time. The latest status is also saved in
`PROGRESS_eivGP.txt` and `PROGRESS_competitors.txt` in the result root;
`_meta/*_status.json` and `logs/` retain machine-readable status and detailed
worker output. Each OI comparator is an independent job, so a slow method does
not hold the other two methods behind it within a fold.

Paper settings initially use 50,000 total iterations per chain, burn 5,000,
HMC O0/O20 and ESS O50. Prediction selects 250 retained states per chain.
`--ess-min=100` optionally continues the saved states until monitored scientific
bulk and tail ESS exceed 100; it does not certify convergence. Explicit `--iter`,
`--burn`, `--stages` and sampler overrides are frozen with the run identity.

### Run one fold and resume

```sh
Rscript tools/oai-eivGP.R --profile=paper --cores=4 --run-folds=1
```

`--run-folds` selects work without changing the frozen fivefold design. Repeat
the command after interruption, change it to `--run-folds=2` for another fold,
or omit it to complete all remaining folds in the same output. Completed
verified jobs are reused. `--folds` instead defines the study scope and must not
be changed inside an existing result root.

With the default paper checkpoint settings, recoverable per-chain states are
saved at 5,000 (end of burn-in). The application then
publishes the first full `fit.rds` at 10,000. The full fit is atomically replaced at
15,000, 20,000, ..., 50,000, retaining all post-burn draws and auxiliary probit S
scores. Each task has `mcmc_progress.csv` and `mcmc_history.csv`; initial per-chain
timing is under its `checkpoints/` directory. The run manifest records the first
checkpoint and checkpoint interval. Continuation follows those recorded settings
and preserves the sampler state, RNG state and retained draws. Exact continuation
requires the same package and numerical environment.

Each result root has a frozen manifest and fitting identity. Reuse and
continuation require matching package modules, numerical settings and input
hashes. A lock records its owner's PID and host. After an interruption, confirm
that the owner and its workers have stopped before clearing a residual lock.

## 3. Pilot and quick checks

`--profile=pilot` uses the same fivefold workflow with a single 1,500-iteration
stage, burn 500 and 250 prediction states per chain. With at least four cores this
is four chains retaining 1,000 states each. Competitor controls remain the regular
paper controls. Pilot and paper profiles use separate output roots, each with
its sampling budget and burn-in recorded in the run manifest.

`--profile=quick` shortens the synthetic check and competitor optimization; it is
not a scientific performance comparison. A local synthetic check exercises the
actual EIV, competitor and report CLIs, including saved-fit reuse:

```sh
Rscript tests/application/check-oai.R --cores=4 --output=results/oai_check
```

No OAI records are used by that check.

On a machine with at least eight available cores, the following additionally
tests two real simultaneous processes (four cores each), numerical agreement
with the sequential check and reuse of saved fits. Use new output directories:

```sh
Rscript tests/application/check-oai-concurrent.R --cores=4 --reference=results/oai_check --output=results/oai_concurrent_check
```

## 4. Saved-only rendering and outputs

```sh
Rscript tools/oai-report.R --profile=paper --cores=4
Rscript tools/oai-eivGP.R --profile=paper --cores=4 --reuse-only=true
Rscript tools/oai-competitors.R --profile=paper --cores=4 --reuse-only=true
```

The combined saved report is `results/oai/oai.html`. Method pages are
`results/oai/reports/eivGP/index.html` and
`results/oai/reports/competitors/index.html`. The editable source is
`paper/applications/oai.Rmd`; `paper/applications/oai.html` is a data-free plan.

| Result block | Contents |
| --- | --- |
| `eiv_foldNN_O*/` | Final fit, predictions and diagnostics; earlier stages stay inside the task |
| `baseline_foldNN/`, `cc_foldNN_O*_METHOD/` | Native OI and complete-case fits and scores |
| `data/`, `_meta/`, `logs/` | Verified private input snapshots, identities, receipts and task logs |
| `reports/` | Method pages, prediction summaries, paper exports and scientific diagnostics |

Use `--output` consistently to select the result root for fitting and reporting.

## 5. Interpretation and provenance

OI predicts all held-out participants from (X,C). CC baselines train only on
revealed calibration U and use (X,U,C); EIV conditions on the same held-out U in
f(X,U). Every CC method uses the identical supplied known-U test subset. The
fixed observation masks are read directly from the prepared protocols.

Training-U recovery uses hidden training rows and may condition on training Y;
prospective recovery uses test C without test Y. O0 has no physical-U error score.
The optional metadata field `imputation_display_fold` selects the main training-U
figure, falling back to the configured representative fold when unavailable.
Response-function and trace displays follow `--representative-fold`.

All finite methods remain in metric rankings. Asterisks identify CORE failures
or exposed flagged optimizer returns; missing diagnostics are not declared
converged. Original results and dataset-specific provenance belong in local
inputs/results, not public source files.
