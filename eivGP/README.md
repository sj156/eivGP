# eivGP 0.4.0.9020

The package fits errors-in-variables Gaussian processes with ordinal inputs.
Public APIs support fitting, prediction, imputation, diagnostics, continuation,
and reporting. The implementation follows the upstream three-module structure:

- `helpers.R`: data preparation, model utilities and analysis workflows, including Ocean execution and core-model calls.
- `mcmc.R`: posterior sampling, fitting and exact continuation.
- `hmc.R`: joint HMC coordinates, warmup adaptation and diagnostics.
- `disk.R`: bounded-memory continuation and immutable draw-history archives.
- `plotting.R`: summaries, tables, figures and report generation, including Ocean reporting and diagnostics.

Package metadata and static-analysis declarations remain in their existing files.
The sampler remains version 0.3.1. Ocean-specific sections have explicit headings;
no separate Ocean implementation scripts are installed.

See the [repository guide](../README.md) for installation,
the [paper guide](../paper/README.md) for the four analysis workflows,
and the [Ocean guide](../paper/README_ocean.md) for all three cohort sizes.
Use `help(package="eivGP")` for function documentation.

Choose `sampler="ESS"` (default), `"HMC"`, or `"HYBRID"`. The hybrid adds
a complete ESS sweep after every 20 HMC iterations; change that frequency
with `sampler_control=list(hybrid_every=20L)`. Continuation preserves the
selected kernel, frequency, frozen adaptation and RNG state.

### Continue without loading the full draw history

```r
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
