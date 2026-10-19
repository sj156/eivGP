\name{plot_mixing}
\alias{plot_mixing}
\title{Visualize saved MCMC diagnostic summaries}
\usage{plot_mixing(object, type = c("overview", "precision", "coverage"),
  rhat_limit = NULL, ess_limit = NULL, mcse_limit = NULL)}
\arguments{
\item{object}{An \code{eivgp_diagnostics} report or data frame with \code{parameter} (or \code{variable}), \code{rhat}, \code{ess_bulk}, \code{ess_tail}, and optionally \code{mcse_sd_ratio}. Optional \code{fit_id}, \code{panel}, and \code{diagnostic_role} columns control grouping. Otherwise replication, scenario and calibration columns are used when present.}
\item{type}{Chain agreement versus bulk ESS, mean precision, or unavailable-value counts.}
\item{rhat_limit, ess_limit, mcse_limit}{Reference lines. NULL uses diagnostic report settings when available, otherwise 1.01, 400, and 0.05. The ESS reference is for bulk ESS only.}}
\value{A customizable ggplot object. Its data contain fit-wise extrema or unavailable counts.}
\description{Plots existing diagnostics without fitting, prediction or integration. Finite extrema are computed separately for each fit and diagnostic role.}
\details{Unavailable diagnostics remain counted in the coverage view, never treated as passing. A point's extrema can come from different parameters. Empty precision information is explicitly labeled. Tail ESS contributes to unavailable counts but is not the overview y axis. These screens do not establish convergence, identification, model adequacy or integration stability.}
\seealso{\code{\link{diagnose_eivgp}}, \code{\link{plot_mcmc}}, \code{\link{save_mixing_plots}}}
\examples{
d <- data.frame(parameter=c("a","b"),rhat=c(1.002,1.03),
  ess_bulk=c(900,150),ess_tail=c(700,NA),mcse_sd_ratio=c(.03,.09))
plot_mixing(d)
plot_mixing(d, "coverage")
}
