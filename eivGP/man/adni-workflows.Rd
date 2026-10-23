\name{adni_workflow_api}
\alias{adni_workflow_api}
\alias{run_adni_application}
\title{Prepared-fold ADNI application workflow}
\description{Package implementation for prepared private ADNI folds, fitting,
scoring, scientific diagnostics and saved-result reports. The application
scripts provide the command-line entrypoints.}
\usage{
adni_workflow_api(repository_dir = NULL)
run_adni_application(repository_dir, params = list(),
  group = c("report", "eivGP", "competitors"))
}
\arguments{
\item{repository_dir}{Repository root. When supplied, the installed source
manifest must match the checkout.}
\item{params}{Named application settings, including profile, input/output paths
and optional numerical overrides.}
\item{group}{Method group or saved-result report.}
}
\value{\code{adni_workflow_api} returns an isolated workflow environment.
\code{run_adni_application} returns the selected workflow result invisibly.}
\details{Participant inputs and empirical reports remain local. Report and
reuse-only actions do not fit missing models or extend chains. This prepared-fold
interface preserves the existing \code{run_adni_analysis} interface.}
