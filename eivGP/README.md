# eivGP 0.4.0.9003

The package fits errors-in-variables Gaussian processes with ordinal inputs.
Public APIs support fitting, prediction, imputation, diagnostics, continuation,
and reporting. The implementation follows the upstream three-module structure:

- `helpers.R`: data preparation, model utilities and analysis workflows, including Ocean execution and core-model calls.
- `mcmc.R`: posterior sampling, fitting and continuation; unchanged from the upstream base.
- `plotting.R`: summaries, tables, figures and report generation, including Ocean reporting and diagnostics.

Package metadata and static-analysis declarations remain in their existing files.
The sampler remains version 0.3.1. Ocean-specific sections have explicit headings;
no separate Ocean implementation scripts are installed.

See the [repository guide](../README.md) for installation,
the [paper guide](../paper/README.md) for the four analysis workflows,
and the [Ocean guide](../paper/applications/README_ocean.md) for all three cohort sizes.
Use `help(package="eivGP")` for function documentation.
