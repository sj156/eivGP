## Install the dependencies needed to install eivGP and render its numerical
## reproduction documents. This deliberately excludes optional published-
## competitor packages; the simulation preflight records their availability.

required <- c(
  "posterior", "TruncatedNormal", # eivGP Imports
  "rmarkdown", "knitr",           # document rendering
  "ggplot2", "dplyr", "tidyr", "patchwork" # tables and figures
)

missing <- required[!vapply(required, requireNamespace, logical(1),
                             quietly = TRUE)]
if (!length(missing)) {
  cat("All required eivGP reproduction dependencies are already installed.\n")
} else {
  cat("Installing required eivGP reproduction dependencies:\n  ",
      paste(missing, collapse = ", "), "\n", sep = "")
  utils::install.packages(missing, repos = "https://cloud.r-project.org")

  still_missing <- required[!vapply(required, requireNamespace, logical(1),
                                     quietly = TRUE)]
  if (length(still_missing)) {
    stop(
      "Installation did not complete for: ", paste(still_missing, collapse = ", "),
      ". Check the R console output and your network/proxy configuration."
    )
  }
  cat("All required eivGP reproduction dependencies are installed.\n")
}
