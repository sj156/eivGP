\name{report_numerical_convergence}
\alias{report_numerical_convergence}
\title{Audit saved numerical convergence and paired accuracy}
\usage{report_numerical_convergence(roots, output)}
\arguments{
\item{roots}{Named character vector of directories containing numerical job registries. Equal names assert verified matching design, data and scoring protocols.}
\item{output}{Separate output directory outside the input roots.}}
\description{Reads saved scores and stage diagnostics without fitting or modifying input artifacts. Reports stage-specific failures, unknown evidence and pairwise accuracy on runs passing all applicable diagnostic screens.}
\value{Invisibly returns a list of status, failure-rate, replicate-score and matched-comparison tables; also writes CSV files.}
\details{Unknown diagnostics never count as passes. Failure denominators and matched sample sizes are explicit. Native adapter errors are distinct from optimizer nonconvergence. Matched results can be affected by selection on successful computation.}
