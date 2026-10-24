\name{continue_study1_results}
\alias{continue_study1_results}
\title{Continue and refresh only deficient saved Study I fits}
\usage{continue_study1_results(path, data_root,
  continuation_dir = file.path(path, "continuation"), ess_min = 100,
  chunk_iter = 5000L, n_cores = 1L, cleanup = TRUE)}
\arguments{
\item{path}{Exact saved Study I run directory.}
\item{data_root}{Frozen Study I data directory.}
\item{continuation_dir}{Directory containing resumable extension states.}
\item{ess_min}{Strict lower bound for scientific-target bulk and tail ESS.}
\item{chunk_iter}{Additional transitions per chain between diagnostic checks.}
\item{n_cores}{Chain worker count. Fits are processed sequentially.}
\item{cleanup}{Remove redundant run checkpoint files after verifying successful results.}}
\description{Resume only fits below the scientific-target ESS requirement, preserve original retained draws, recompute their prediction, mean, surface and input scores, and refresh replication bundles and combined exports. Other fits and their scores are reused.}
\value{Updated run directory and continuation audit.}
\details{Frozen-data checksums, target panels and retained draw prefixes are checked before replacement. Reports can subsequently be rendered with report_study1_eivgp(). Final fit terminal states are retained.}
