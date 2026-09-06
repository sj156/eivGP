args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
rmd_path <- normalizePath(
  file.path(dirname(script_path), "create-eivmixgp.Rmd"),
  mustWork = TRUE
)
## litr 1.2.0's HTML hyperlink post-processor is not vector-safe when a
## package exposes a few hundred internal functions.  Markdown still exercises
## the complete literate build, tests, roxygen pass, versioning, and hashing,
## while avoiding that presentation-only post-processing failure.
litr::render(
  rmd_path,
  output_format = rmarkdown::md_document()
)

## CRAN reserves arbitrary DESCRIPTION extensions under Config/*.  Preserve
## litr's provenance values there instead of shipping its legacy top-level
## fields, which R CMD check reports as unknown.
package_dir <- normalizePath(
  file.path(dirname(rmd_path), "..", "package-build", "eivmixgp"),
  mustWork = TRUE
)
description_path <- file.path(package_dir, "DESCRIPTION")
description <- read.dcf(description_path)
description_fields <- as.list(description[1L, , drop = TRUE])
description_fields[["Config/litr/version"]] <-
  description_fields[["LitrVersionUsed"]]
description_fields[["Config/litr/id"]] <- description_fields[["LitrId"]]
description_fields[c("LitrVersionUsed", "LitrId")] <- NULL
write.dcf(
  as.data.frame(
    description_fields,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  file = description_path,
  width = 80
)
