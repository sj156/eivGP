# Bounded mixing case study

From the repository root:

```sh
Rscript experiments/case_study_mixing.R --smoke
Rscript experiments/case_study_mixing.R 2500 4
```

The arguments are additional iterations per chain, worker cores, and an optional
new output directory. Default: 2,500 additional transitions per chain, four
workers, a timestamped directory under `reproduction/case-studies/`.
The script does not install packages. Use the repository dependencies already
required for fitting and diagnostics. `EIVGP_REPO` optionally sets the repo path.

This case is Study I eta1_balanced, replication 1, 20 calibrated observations,
saved in development run `study1-development-ef651d075f78`. The script requires
that exact saved fit. It is a repository-specific legacy checkpoint adapter,
not a new package continuation API. Public fits should use `continue_eivgp()`.

The saved fit has 1,250 retained draws per chain. The default adds 2,500 for
3,750 retained draws per chain (15,000 total). No new burn-in, thinning, prior
changes, or block-size changes are introduced. Terminal states and RNG streams
are resumed. Original prefixes, data, priors, and the input-file hash are checked.
The smoke test verifies that a two-step continuation equals two one-step calls.

Outputs include the extended fit RDS, before/after diagnostic tables, a joined
comparison CSV, per-chain moments/quantiles, dictionary occupancies, and trace
PDFs. Unvisited dictionary states and undefined ESS/R-hat remain visible.
The complete extended fit is saved before after-run target diagnostics so a
reporting failure does not erase the completed sampling work.

Both reports use the same five calibration-anchored training locations for
`m(x,c)` and `f(x,u)`, all retained draws, 32 latent integration points, and
the same integration seed. These are a diagnostic panel, not a held-out accuracy
experiment or whole-domain certificate. Their finite integration error is not
included in MCMC MCSE. Later sensitivity checks should increase integration
points or change the integration seed.

Compare chain-specific means and traces, especially signal fraction, noise,
cutoffs, and dictionary occupancy. Check whether scientific-target R-hat moves
toward one and mean-specific MCSE decreases. Increased ESS alone does not prove
that modes have been explored. No diagnostic gate blocks completion and no
additional run is triggered automatically. Saved before/after summaries use
their existing timing metadata; `extension_time.rds` records this extension's
wall time separately, so do not interpret after-fit timing as total lifetime cost.
