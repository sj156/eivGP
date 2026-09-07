# Literate eivGP package build

`create-eivGP.Rmd` is the sole active build source. It copies canonical model
modules from `../codes` into `../eivGP/R`, then creates package metadata,
documentation and tests. `../eivGP` is the only current installable package.

Run from the repository root:

```sh
Rscript --vanilla litr/render-package.R
R CMD build eivGP
R CMD check --no-manual eivGP_0.3.1.tar.gz
```

The wrapper resolves absolute paths, so the output location does not depend
on the caller's working directory. Model fixes belong in `codes/`; generated
package modules must be byte-identical to their canonical files. The build
sets their load order explicitly so the current sampler definitions are active.

Study generators and orchestration are repository code. The package does not
install `00_synthetic_data.R`, `00_experiment_runner.R`, or `inst/experiments`.
Current target/workflow tests and reusable utility tests are installed with the
package. Retired literate sources and their obsolete tests have been removed
from the working tree; previously committed versions remain in Git history.
