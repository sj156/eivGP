\name{continue_eivgp_until_ess}
\alias{continue_eivgp_until_ess}
\title{Continue posterior chains to a scientific-target ESS requirement}
\usage{continue_eivgp_until_ess(object, X, C, U = NULL, ess_min = 100,
  chunk_iter = 5000L, max_extensions = Inf, n_latent = 64L,
  seed = 481517L, parallel = FALSE, n_cores = 1L, checkpoint_file = NULL,
  adaptive = FALSE)}
\arguments{
\item{object}{An eivgp_fit or supported internal Study I or Study II fit.}
\item{X, C, U}{Fixed diagnostic locations in original units; U is optional.}
\item{ess_min}{Strict lower bound for both bulk and tail ESS.}
\item{adaptive}{Estimate each additional budget from the weakest ESS, bounded between chunk_iter and four times chunk_iter; reassess after each extension.}

\item{chunk_iter}{Additional transitions per chain between checks.}
\item{max_extensions}{Optional operational limit, otherwise continue until success.}
\item{n_latent, seed}{Fixed common-random-number integration settings.}
\item{parallel, n_cores}{Chain execution controls.}
\item{checkpoint_file}{Optional atomic RDS destination after every extension and check.}}
\value{A list containing fit, target_series, diagnostics, history, settings and success.}
\description{Continue existing chains without dropping retained draws or repeating warmup until all monitored scientific targets have bulk and tail ESS strictly above the threshold. At least four chains and full retained windows are required.}
\details{Undefined ESS never passes. R-hat and MCSE remain reported, but do not control this precision rule. Success does not establish convergence; inspect traces. An operational limit returns success=FALSE when unmet. The checkpoint retains complete draws and terminal states even if diagnostics are interrupted. Read it with read_fit() and pass its fit back to resume.}
