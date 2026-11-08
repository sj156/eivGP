\name{continue_eivgp_disk}
\alias{continue_eivgp_disk}
\title{Continue MCMC with complete draw history on disk}
\description{Advance the existing chains without loading or joining the retained draw history during sampling.}
\usage{
continue_eivgp_disk(path, n_iter, chunk_iter = 1000L, parallel = TRUE,
  n_cores = NULL, verbose = FALSE, progress = NULL)
}
\arguments{
\item{path}{Existing public fit checkpoint path.}
\item{n_iter}{Additional transitions per chain; zero imports an existing fit only.}
\item{chunk_iter}{Maximum new retained draws per chain buffered in memory.}
\item{parallel, n_cores}{Chain execution controls.}
\item{verbose}{Print committed progress.}
\item{progress}{Optional callback receiving the small committed index.}
}
\value{Invisibly, a small disk-backed fit index. The checkpoint is updated atomically.}
\details{The original fit is preserved at \code{path.draws/legacy-fit.rds}. Import reads one legacy fit once; subsequent sampling uses terminal Markov states, RNG states, unchanged tuning, absolute iteration numbers and bounded draw buffers. Immutable draw blocks and HMC diagnostic traces remain on disk. No thinning or new warmup is performed. The archive index commits after all chains' blocks are durable; uncommitted blocks are recomputed after interruption. The application executor must provide a single-writer lock. Keep the checkpoint and its adjacent \code{.draws} directory together when moving results. \code{read_fit(path)} explicitly reconstructs all history for compatible analysis, requiring memory proportional to that history. Its optional \code{draw_ids} and \code{fields} arguments allow selective reads; only auxiliary S fields may be omitted from diagnostic fits. Selected fits cannot be continued with the ordinary in-memory API; continue from the archive path instead.}
\seealso{\code{\link{continue_eivgp}}, \code{\link{read_fit}}}
