#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("remotes", "knitr", "jsonlite", "ggplot2", "patchwork", "dplyr", "tidyr", "posterior", "kergp", "LVGP", "EzGP")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
commit <- "a3aa0f697f97d4bdf027507ac3c9b952644affb9"
pkg <- if (requireNamespace("eivGP", quietly = TRUE)) utils::packageDescription("eivGP") else NULL
if (is.null(pkg) || !identical(pkg$Version, "0.3.1") || !identical(pkg$RemoteSha, commit))
  remotes::install_github("sj156/eivGP", subdir = "eivGP", ref = commit,
                          upgrade = "never", dependencies = NA, force = TRUE)
writeLines(capture.output(sessionInfo()), "setup-sessionInfo.txt")
