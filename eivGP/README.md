# eivGP 0.4.0.9023

The package fits errors-in-variables Gaussian processes with ordinal inputs.
Public APIs support fitting, prediction, imputation, diagnostics, continuation,
and reporting. The implementation is organized by responsibility:

- `helpers.R`: data preparation, shared model utilities and numerical-study workflows.
- `mcmc.R`: posterior sampling, fitting and exact continuation.
- `hmc.R`: joint HMC coordinates, warmup adaptation and diagnostics.
- `disk.R`: bounded-memory fresh fitting, continuation and immutable draw-history archives.
- `plotting.R`: summaries, tables, figures and shared report/plot rendering.

Package metadata and static-analysis declarations remain in their existing files.
The sampler remains version 0.3.1. Application-specific orchestration stays in
the flat `tools/APP-tools.R` files; the package supplies reusable calculations.
Unused pre-0.3.1 sampler and pilot-baseline helpers have been removed. Public
APIs, dynamically selected study constructors/reports and checkpoint readers
remain available.

See the [repository guide](../README.md) for installation,
the [paper guide](../paper/README.md) for the four analysis workflows,
and the [Ocean guide](../paper/README_ocean.md) for all three cohort sizes.
Use `help(package="eivGP")` for function documentation.

Choose `sampler="ESS"` (default), `"HMC"`, or `"HYBRID"`. The hybrid adds
a complete ESS sweep after every 20 HMC iterations; change that frequency
with `sampler_control=list(hybrid_every=20L)`. Continuation preserves the
selected kernel, frequency, frozen adaptation and RNG state.

### Fit and continue with draws on disk

```r
index <- fit_eivgp(X, y, C, U_obs = U_obs, engine = "multivariate",
  latent_dim = 1L, ident = "none", n_iter = 50000L, burn = 5000L,
  n_chains = 4L, n_cores = 4L,
  storage = "disk", path = "fit.rds", chunk_iter = 1000L)

# Import an existing fit once, then add 10,000 transitions per chain.
continue_eivgp_disk("fit.rds", n_iter = 10000L,
                    chunk_iter = 1000L, parallel = TRUE, n_cores = 4L)

# Analysis is explicit; this call reconstructs the complete history in RAM.
fit <- read_fit("fit.rds")
```

The checkpoint becomes a small terminal-state index. Its adjacent
`fit.rds.draws/` directory holds every retained draw, the complete HMC
diagnostic trace, checksums and the preserved original checkpoint. Keep both
together when moving results. Initial import reads one legacy fit; subsequent
sampling keeps only the current states and at most `chunk_iter` new draws per
chain in RAM. Interrupted writes resume from the last committed block with the
same RNG states and absolute iteration numbers. No warmup or thinning is added.

Fresh disk fits run the unchanged sampler through its original warmup and at
most one retained block, then use disk-backed continuation. Configured warmup
recovery checkpoints remain supported. The initial HMC warmup trace is held
until that first publication; retained-history memory is bounded by the block
size rather than the total iteration budget. Existing checkpoint paths are
refused; resume them with `continue_eivgp_disk()`.

ADNI, OAI and Ocean real-data workflows use disk storage by default. Their
sampling phases do not materialize past draws. Full-history evaluation and
scientific ESS checks run separately. ADNI and Ocean use one resident fit;
OAI allows one evaluation per four cores, capped at eight jobs. Each evaluation
still needs memory proportional to its parameter history. Automatic scientific
ESS uses the same check and adaptive extension policy, storing its target
series on disk and releasing the analysis fit before further sampling.
Direct R calls retain `storage="memory"` as the compatibility default.

`read_fit(path, draw_ids=...)` can reconstruct selected draws (at least two per
chain). `read_fit(path, fields=...)` can omit auxiliary `S` draws while retaining
complete parameter histories for diagnostics. These selected analysis fits are
not continuation checkpoints; continue from the archive path. A full-history
read has memory cost proportional to the stored history. Application executors
provide the single-writer lock required for disk continuation.

### Dependency scope

The numerical API loads Rfast, Rcpp, posterior and TruncatedNormal (plus standard R libraries);
RcppArmadillo is needed to build the native backend. Plotting and report packages
are optional and loaded only when those operations are requested. Small internal
helpers implement ordinary data-frame row binding, long/wide metric reshaping,
and calibration-grid expansion. They use explicit column names and unique
wide-table keys, rather than providing general replacements for a data toolkit.
Scientific diagnostics, minimax-tilting sampling, native linear algebra, JSON
parsing and provenance hashes retain their established package implementations.
