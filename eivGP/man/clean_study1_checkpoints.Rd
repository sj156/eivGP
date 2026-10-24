\name{clean_study1_checkpoints}
\alias{clean_study1_checkpoints}
\title{Remove redundant checkpoint files after Study I success}
\usage{clean_study1_checkpoints(path, ess_min = 100)}
\arguments{\item{path}{Saved Study I run directory.}\item{ess_min}{Strict scientific-target ESS lower bound.}}
\description{Require complete scores and successful scientific-target diagnostics for every expected fit. Before deleting checkpoint files, verify that final fit histories and terminal states exist. Active locks prevent cleanup.}
\value{File removal audit, also saved under config/checkpoint_cleanup.csv.}
\details{Only files under the run's checkpoints directory are deleted. Final fit artifacts remain available for diagnostics and continuation. Failure leaves checkpoints intact.}
