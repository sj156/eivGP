if (!requireNamespace("eivGP", quietly = TRUE) ||
    as.character(utils::packageVersion("eivGP")) != "0.3.1")
  stop("Both applications require eivGP version 0.3.1. Run codes/setup_real_data.R.")
