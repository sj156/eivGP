# Manuscript follow-up simulations

`study1.R` and `study2.R` are complete fresh-run entrypoints. Our additions use
`temp-study1-followup.R` and `study2-eivGP.R` and `study2-competitors.R`, which call the same package
implementation and can be deleted after the follow-up work is complete.

The follow-up grid does not rerun either original core experiment or the
original categorical-only competitors. It adds 600 Study I complete-case fits
(both regimes, calibration 20/50) and 400 Study II fits (100 CUT plus 300
complete-case competitors), and evaluates 100 saved Study II FULL fits at
matched independent test locations. Original study results remain unchanged.

## Study II paper replication

Use [the separated Study II workflow](../study2-README.md) for current commands,
remaining-job selections, and reporting. The old combined temporary Study II
launcher is retired. No pilot or integration variant is required to regenerate
the existing paper displays. FULL is saved-fit evaluation, CUT is the secondary
feedback comparison, and the competitor lane fits the published complete-case
adaptations independently. Keep their output roots separate across machines.

Study I continues to use `temp-study1-followup.R` with the documented original
budgets; its workflow is unchanged by the Study II split.

### Automatic replication parallelism

For the follow-up entrypoints, `--cores=N` is the total CPU budget. On Mac and
Linux, the runner dynamically schedules independent replications with up to four
chain processes per replication: `--cores=56` runs up to 14 replications at once.
Methods and calibration levels within each replication stay sequential so they
can safely reuse measurement fits. Study I runs its two regimes sequentially.
Prediction and optimization phases may use fewer cores than MCMC; utilization is
not guaranteed to stay at 100%. More simultaneous replications require more RAM.

Keep the same `--reps`, `--output`, and other scientific options when restarting.
The scheduler retains the original batch paths and skips completed scores. It
accepts the explicitly recorded scheduling-only predecessor hashes, while still
rejecting changes to inference code, dependencies, data, or scientific settings.
Stop the old run before updating/restarting; never start a second writer in its
output directory. Abrupt termination may leave locks: confirm all old parent and
child processes have exited before removing a stale lock.

Use single-threaded numerical libraries when launching R to avoid nested BLAS
oversubscription:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

No external `xargs` pool is needed. Install the matching updated GitHub package
before restarting these entrypoints. Existing processes retain their old scheduler.
