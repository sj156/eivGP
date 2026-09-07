args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
rmd_path <- normalizePath(
  file.path(dirname(script_path), "create-eivGP.Rmd"),
  mustWork = TRUE
)
## R CMD check reserves custom DESCRIPTION fields under Config/*. Configure
## litr before rendering so its hash covers the final DESCRIPTION, rather
## than changing metadata after hashing and breaking subsequent builds.
## This wrapper is already an isolated Rscript process; keeping the render in
## this session ensures the metadata-field adapters apply to setup and cleanup.
utils::assignInNamespace("description_litr_hash_field_name",
  function() "Config/litr/id", ns = "litr")
utils::assignInNamespace("description_litr_version_field_name",
  function() "Config/litr/version", ns = "litr")

## Markdown exercises the complete literate build and documentation while
## avoiding the HTML post-processor's large-namespace hyperlink limitations.
litr::render(
  rmd_path,
  fresh_session = FALSE,
  output_format = rmarkdown::md_document()
)
