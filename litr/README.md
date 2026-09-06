# Literate package build

`create-eivmixgp.Rmd` is a transitional `litr` 0.9.x build document. The
canonical code remains in `../codes`; rendering copies the audited modules into
`../package-build/eivmixgp/R` and generates package metadata, `README.md`,
`NEWS.md`, reference documentation, and tests.
The manuscript experiment scripts are installed under `inst/experiments`, so
`run_mixedgp_experiment()` provides the same Study I/Study II launcher from a
development checkout or an installed package.

Run from this directory:

```sh
Rscript render-package.R
```

The wrapper supplies an absolute input path to `litr::render()`, which keeps
the generated package location stable regardless of the caller's working
directory.

Do not patch files under `../package-build/eivmixgp` directly. Change the
canonical modules or `create-eivmixgp.Rmd`, render, and verify that the copied
modules are byte-identical to their canonical sources.

Remaining release follow-up:

1. split plotting/reporting helpers from the computational namespace;
2. move a small demonstration dataset into package data and retain the full
   publication datasets under `data-synthetic/` with their checksum manifest;
3. run `devtools::check(document = FALSE)` on macOS, Linux, and Windows.
