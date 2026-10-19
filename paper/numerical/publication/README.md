# Manuscript follow-up simulations

`study1.R` and `study2.R` are complete fresh-run entrypoints. Our additions use
`temp-study1-followup.R` and `temp-study2-followup.R`, which call the same package
implementation and can be deleted after the follow-up work is complete.

The follow-up grid does not rerun either original core experiment or the
original categorical-only competitors. It adds 600 Study I complete-case fits
(both regimes, calibration 20/50) and 400 Study II fits (100 CUT plus 300
complete-case competitors), and evaluates 100 saved Study II FULL fits at
matched independent test locations. Original study results remain unchanged.

## Two-machine allocation

From the repository root, first install this source on each host:

```sh
Rscript tools/setup_numerical.R
```

Set machine-local directories explicitly. Mac mini uses replications 1:25;
Linux uses 26:50. Each machine needs the two saved Study II FULL archives per
replication at calibration 20 and 50 (50 archives per host). No previous
Study I artifacts, original competitor archives, or old pilot folders are
required to execute the new jobs. Frozen synthetic data are regenerated and
FULL matching checks verify data, calibration, priors and model.

```sh
# Mac mini (run from the checkout root)
Rscript paper/numerical/temp-study1-followup.R --action=run --profile=paper --cores=4 --reps=1:25 --output="$HOME/eivGP-followup-20261019/mac-mini/study1" --data-root="$HOME/eivGP-followup-20261019/mac-mini/data/study1"
Rscript paper/numerical/temp-study2-followup.R --action=run --profile=paper --cores=4 --reps=1:25 --output="$HOME/eivGP-followup-20261019/mac-mini/study2" --data-root="$HOME/eivGP-followup-20261019/mac-mini/data/study2" --full-fit-root="$HOME/eivGP-followup-20261019/mac-mini/full-fits"

# Linux (run from the checkout root)
Rscript paper/numerical/temp-study1-followup.R --action=run --profile=paper --cores=4 --reps=26:50 --output="$HOME/eivGP-followup-20261019/linux/study1" --data-root="$HOME/eivGP-followup-20261019/linux/data/study1"
Rscript paper/numerical/temp-study2-followup.R --action=run --profile=paper --cores=4 --reps=26:50 --output="$HOME/eivGP-followup-20261019/linux/study2" --data-root="$HOME/eivGP-followup-20261019/linux/data/study2" --full-fit-root="$HOME/eivGP-followup-20261019/linux/full-fits"
```

Replace the home-directory output paths with a large non-synced disk if needed.
Run one study at a time per machine initially. Use `--action=plan` to inspect
without fitting. Repeat identical commands to resume compatible completed
work; do not change batches, code or budgets within an existing output root.
Do not remove stale locks until their owner is confirmed stopped. Old pilot
caches are not silently imported across incompatible protocols.

The current budgets remain provisional: Study I 12,000/4,000 total/warmup,
16 outer states, 128 latent and 64 response draws; Study II measurement
100,000/20,000, 32 outer states, conditional GP 4,000/1,000 with four selected
states per chain, 128 latent, 16 GP-response and 64 competitor-response draws.
All use four chains. Completing this grid supplies the missing comparisons;
it does not certify posterior uncertainty or integration adequacy. Targeted
integration checks are available through `--stage=integration --reps=1,25`,
with separate outputs. These variants can require new conditional fits, so do
not launch them as an automatic repeat of the complete grid. Preserve warnings,
failed fits, diagnostics and timing. Replication MCSE omits within-fit error.

For final reporting, combine the new scores with existing core results by
study/regime/replication/calibration/method; reject duplicates rather than
averaging host-level means. Integration variants are not new replications.

Exploratory sampler experiments are kept outside the repository; they are
not prerequisites for these follow-up commands. Superseded shell launchers
have been retired outside the repository. No application scripts are changed.
