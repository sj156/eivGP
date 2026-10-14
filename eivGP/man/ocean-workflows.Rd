\name{ocean-workflows}
\alias{run_ocean_analysis}
\alias{ocean_settings}
\alias{ocean_workflow_api}
\title{Ocean application workflows and report definitions}
\description{Run or report frozen 200- or 500-profile Ocean designs using existing model and competitor implementations.}
\usage{
run_ocean_analysis(data_dir, output_dir, action = "plan", profile = "paper",
  cores = 4L, repository_dir = NULL, settings = NULL, import_dirs = character(),
  cc_output_dir = "", diagnostics_dir = "", representative_fold = 1L,
  diagnostics_workers = 1L, ...)
ocean_settings(profile = "paper", action = "plan", folds = 1:5,
  cells = c("O0", "O20", "O50"),
  modules = c("eiv", "original_baselines", "continuation",
    "complete_case_baselines", "prediction", "diagnostics", "surface_comparison"),
  task_ids = character(), cores = 4L, run_id = "local")
ocean_workflow_api()
}
\arguments{
\item{data_dir}{Frozen cohort, fold and mask input directory.}
\item{output_dir}{Separate output root; new identities require a new root.}
\item{action}{Plan, run or report; report never refits.}
\item{profile}{Paper, complete short-budget pilot, or single-fold quick subset.}
\item{cores, diagnostics_workers}{CPU budgets.}
\item{repository_dir}{Repository root, defaulting to the working directory.}
\item{settings}{Optional resolved settings.}
\item{import_dirs}{Compatible task roots to import.}
\item{cc_output_dir, diagnostics_dir}{Optional CC and diagnostic destinations.}
\item{representative_fold}{Main-figure fold selected in advance.}
\item{folds, cells, modules, task_ids}{Task selectors.}
\item{run_id}{Log identifier.}
\item{...}{Additional ocean_settings arguments.}
}
\value{Settings or analysis result lists. ocean_workflow_api returns package definitions for the notebook.}
\details{OI, CC and MI protocols remain separate. Frozen masks and data determine sample counts; no hidden input is inserted into CC prediction. Numerical flags are retained. The MCMC engine, native adapters, model and priors are unchanged.}
