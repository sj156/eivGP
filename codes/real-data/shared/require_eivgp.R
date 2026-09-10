if (!requireNamespace("eivGP", quietly = TRUE) ||
    as.character(utils::packageVersion("eivGP")) != "0.3.1")
  stop("Both applications require eivGP version 0.3.1. Run applications/setup.R.")
