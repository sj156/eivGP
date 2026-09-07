# eivmixgp 0.2.1 — local release

Historical release notes only. The current package is `eivGP` 0.3.0; use
`README.md` and `litr/render-package.R`. Paths and commands below record the
retired release and are not current installation instructions.

Canonical source: `codes/`; literate build: `litr/create-eivmixgp.Rmd`;
installable package: `package-build/eivmixgp/`.

Implemented:

- Exported `eivgp_run_settings(core_budget, pending_datasets)`.
- Four chains with 500 warmup + 1,250 retained transitions per chain.
- Core-budgeted dataset and chain concurrency in both numerical drivers.
- Explicit continuation preserved; no automatic extensions.
- R Markdown launcher accepts `MIXEDGP_CORE_BUDGET` (12 or 16, for example).

Validation on 2026-09-06:

- Literate build and source tarball build succeeded.
- `R CMD check --no-manual`: Status OK; 302 passing test expectations,
  zero failures, warnings, or skips in the test suite.
- Both development plan documents rendered successfully.
- Driver allocation validation and real nested-fork PID checks passed.
- No full numerical or publication experiments were launched.

This is not completion of the full development orchestration brief. See
`DEVELOPMENT.md` for remaining work: runtime-failure recovery, time caps,
separate fit/evaluation caches, and the expanded design.
The follow-up warning-policy revision removes diagnostic gate stops in both
modes, keeps flagged fit bundles and guidance, and allows Study II PI-GP/CC-GP
evaluation after measurement diagnostic warnings. Missing competitors are
nonfatal by default. Publication iteration budgets remain longer.
Tests execute the actual diagnostic branches with deliberately failed checks
and verify warning delivery and fit preservation; invalid configurations still
raise errors. Cache schemas changed to avoid stale suppressed-result reuse.

Review the local Git diff before committing and pushing. Install locally with:

```sh
R CMD INSTALL package-build/eivmixgp
Rscript --vanilla -e 'stopifnot(packageVersion("eivmixgp") == "0.2.1")'
```
