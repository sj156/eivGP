\name{run_paper_study}
\alias{run_paper_study}
\title{Run a complete numerical paper study}
\usage{run_paper_study(study = c("study1", "study2"), args = list(), repo_root = getwd())}
\arguments{
\item{study}{Univariate or multivariate study.}
\item{args}{Named CLI-style character options. Component defaults to all; core runs the original grid and mechanism runs the targeted Study II comparisons. Action defaults to plan.}
\item{repo_root}{Checkout containing renv.lock and numerical scripts.}}
\description{Orchestrates frozen data, fitting, diagnostics and scoring through the permanent study entrypoints. Study II includes complete-case and matched CUT comparisons. Reports retain pending diagnostics and do not establish publication readiness. See paper/numerical/README.md for options and retention policy.}
\value{Invisibly, component results.}
