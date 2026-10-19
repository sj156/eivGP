\name{plot_mcmc}
\alias{plot_mcmc}
\title{Plot retained chains with ggplot2}
\usage{plot_mcmc(object, parameters = NULL,
  type = c("trace", "rank", "acf", "running_mean"),
  max_points = 2000L, lag_max = 100L, bins = 20L)}
\arguments{
\item{object}{An eivGP fit with chain-resolved draws, or a uniquely named list of functionals. Each functional is a list of at least two numeric chain vectors with at least four finite retained draws each. Named functionals can include previously computed scientific targets.}
\item{parameters}{Names to display. NULL selects the first six available raw functionals; this is not a worst-case selection.}
\item{type}{Trace, pooled-rank balance, autocorrelation or cumulative mean.}
\item{max_points}{Approximate maximum displayed points per chain (plus final endpoint). Display subsampling does not change calculations or saved draws.}
\item{lag_max}{Maximum ACF lag, reduced for shorter chains.}
\item{bins}{Number of pooled-rank bins.}}
\value{A customizable ggplot object with plotting data.}
\description{Draws retained-chain diagnostics without refitting or computing new predictions or target integrations.}
\details{Running means, ranks and ACF use all supplied retained draws. Chain lengths may differ; rank frequencies are normalized within chain. Ties receive average ranks, so discrete/constant rank plots require care. Constant-chain ACF is undefined and explicitly noted. No draws with missing/nonfinite values are silently removed. Indexing starts at the first retained draw, not the original sampler iteration. Warmup must not be included in supplied vectors. Plotting two chains is allowed but does not satisfy the package's four-chain diagnostic screen.}
\seealso{\code{\link{plot_mixing}}, \code{\link{diagnose_eivgp}}}
\examples{
set.seed(1)
s <- list(a=lapply(1:4,function(i) rnorm(200)))
plot_mcmc(s, "a")
plot_mcmc(s, "a", type="rank")
}
