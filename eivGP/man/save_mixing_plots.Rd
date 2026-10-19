\name{save_mixing_plots}
\alias{save_mixing_plots}
\title{Save default mixing figures from a results directory}
\usage{save_mixing_plots(path, output_dir = file.path(path, "mixing-plots"),
  formats = c("pdf", "png"))}
\arguments{
\item{path}{Directory containing saved study diagnostic CSV files. Searched recursively.}
\item{output_dir}{Separate destination for generated figures and a source manifest. Existing figures with the same names are replaced.}
\item{formats}{One or both of pdf and png.}}
\value{Invisibly, a character vector of figure paths. An empty vector with a warning if parameter-level tables are unavailable.}
\description{Saves overview, precision and coverage plots using \code{plot_mixing}. Called by both paper study workflows after report or run unless disabled.}
\details{Reads per-parameter and target MCMC CSVs, measurement diagnostics, conditional GP diagnostics and FULL diagnostics. Combined diagnostic-detail CSVs are used only when separate tables are absent from their directory, preventing duplicate counts. AppleDouble metadata files are ignored. All discovered saved runs are kept distinct by their relative directory; point selection is not restricted to passing fits. A CSV manifest records source paths and MD5 hashes. Historical run-level gate summaries cannot substitute for parameter-level diagnostics. No fitting, continuation or integration occurs; large chain artifacts are not loaded.}
\seealso{\code{\link{plot_mixing}}, \code{\link{plot_mcmc}}, \code{\link{run_paper_study}}}
