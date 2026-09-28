# No fitting: reader profiles, cache safeguards, and public data integrity.
source("codes/simulation_helpers.R")
source("codes/simulations/notebook_workflow.R")
source("codes/reporting/notebook_helpers.R")
source("codes/applications/ocean_silicate/workflow.R")
repo <- normalizePath(".")
notebooks <- list.files("replication", pattern = "[.]Rmd$", recursive = TRUE)
stopifnot(identical(sort(notebooks), c("01_numerical_experiments.Rmd", "02_ocean.Rmd")))
for (study in c("study1", "study2")) {
  paper <- replication_config(study, "paper", 4L, repo)
  quick <- replication_config(study, "quick", 2L, repo)
  stopifnot(paper$mcmc$n_iter == 20000L, paper$mcmc$burn == 5000L,
    paper$mcmc$n_chains == 4L, all(vapply(paper$cells, `[[`, 1L, "n_rep") == 50L),
    all(vapply(paper$cells, `[[`, 1L, "n") == 100L),
    all(vapply(paper$cells, `[[`, 1L, "n_test") == 100L),
    quick$mcmc$n_iter == 120L, quick$mcmc$n_chains == 1L,
    !identical(paper$data_root, quick$data_root), !identical(paper$output_root, quick$output_root))
}
cohort <- ocean_read_inputs("data/ocean_silicate")
stopifnot(nrow(cohort) == 200L, identical(as.integer(table(cohort$oxygen_class)), rep(50L, 4L)))
tmp <- tempfile("notebook-check-"); dir.create(tmp)
identity <- list(profile = "quick", data = "fixture")
notebook_cache_identity(file.path(tmp, "run"), identity)
notebook_cache_identity(file.path(tmp, "run"), identity)
stopifnot(inherits(try(notebook_cache_identity(file.path(tmp, "run"), list(profile = "paper")), silent = TRUE), "try-error"))
depth <- sink.number()
answer <- notebook_logged(file.path(tmp, "execution.log"), {message("message recorded");warning("warning recorded");42L})
stopifnot(answer == 42L, sink.number() == depth,
  identical(attr(answer, "notebook_warnings"), "warning recorded"),
  any(grepl("message recorded", readLines(file.path(tmp, "execution.log")))))
try(notebook_logged(file.path(tmp, "failure.log"), stop("expected failure")), silent = TRUE)
stopifnot(sink.number() == depth)
message("Two notebooks, paper/quick profiles, public data, logging, and cache safeguards passed.")
