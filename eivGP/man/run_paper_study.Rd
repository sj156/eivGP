\name{run_paper_study}
\alias{run_paper_study}
\title{Run a complete numerical paper study}
\usage{run_paper_study(study = c("study1", "study2"), args = list(), repo_root = getwd())}
\arguments{
\item{study}{Univariate or multivariate study.}
\item{args}{Named CLI-style character options. Component defaults to all; core runs the original grid, complete-case runs the Study I additions, and mechanism runs the targeted Study II comparisons. Action defaults to plan. The mixing-plots option defaults to true and saves three figures in each selected component's mixing-plots folder after report or run; set it to false to skip plotting. No diagnostic plotting occurs during plan, data or fit. The gp-response-latent and cc-response-latent options separately control prospective input integration for response scores (defaults 1 and 4); latent controls observable-mean integration. The pilot profile selects two replications and two calibration levels for Study II mechanism, or one replication in both regimes at calibration 50 for Study I complete-case.}
\item{repo_root}{Checkout containing renv.lock and numerical scripts.}}
\description{Orchestrates frozen data, fitting, diagnostics and scoring through the permanent study entrypoints. Both studies include complete-case published predictors; Study II also includes matched CUT comparisons. Reports retain pending diagnostics and do not establish publication readiness. See paper/numerical/README.md for options and retention policy.}
\value{Invisibly, component results, with generated figure paths in the mixing_plots attribute when plotting is enabled.}

\section{Publication completion}{Use component = "publication" with stage = "pilot", "integration", or "full" in args. The default action is "plan". Both studies cover calibration sizes 20 and 50; Study I covers both regimes. Study II run requires a full-fit-map for saved FULL evaluation. Settings are provisional; stages never automatically advance. See paper/numerical/publication/README.md in the repository.}

\section{Temporary follow-up entrypoints}{The repository temp-study1-followup.R and temp-study2-followup.R scripts select component = "follow-up", which defaults to the complete missing-comparison grid without rerunning the core studies. Study II accepts full-fit-root to discover the saved FULL archives, or an explicit full-fit-map. Set reps to disjoint batches on separate hosts and keep output roots separate.}
