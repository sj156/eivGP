\name{report_study1_eivgp}
\alias{report_study1_eivgp}
\title{Report saved Study I EIV-GP results}
\usage{report_study1_eivgp(path, output_dir = NULL, data_root = NULL,
  ess_min = 100, mixing_plots = TRUE)}
\arguments{
\item{path}{Exact saved run, or a parent containing exactly one Study I run.}
\item{output_dir}{Separate report output directory. Defaults to a sibling study1-report directory.}
\item{data_root}{Frozen Study I input root, required to reconstruct historical traces.}
\item{ess_min}{Scientific-target bulk and tail ESS must both be strictly greater than this value.}
\item{mixing_plots}{Generate three small sets of chain traces, serially.}
}
\description{Reads saved scores and full-window scientific-target diagnostics, writes report tables and a report.rds, and optionally plots selected chains. Never fits or continues MCMC. Historical target traces are reconstructed only for selected retained fits, using verified frozen inputs and original diagnostic settings.}
\value{A list with score summaries, scientific-target diagnostics, trace paths, completeness and provenance.}
