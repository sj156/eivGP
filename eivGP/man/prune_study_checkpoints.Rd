\name{prune_study_checkpoints}
\alias{prune_study_checkpoints}
\title{Verify and Prune Redundant Study II Conditional GP Checkpoints}
\usage{
prune_study_checkpoints(output_dir, dry_run = TRUE, audit_dir,
  progress_every = 25L)
}
\arguments{
\item{output_dir}{Study output directory containing \code{checkpoints/} and \code{mechanism/}.}
\item{dry_run}{Logical. The default audits without modifying the output tree.}
\item{audit_dir}{Required fresh local directory outside the output tree for a CSV manifest and per-fit RDS receipts.}
\item{progress_every}{Positive integer number of directories between progress messages.}
}
\description{
Verifies legacy Study II CUT restart files against published complete conditional fits.
Candidate filenames and seeds are used only for lookup. All model inputs, settings,
chain draws, initial and terminal states, RNG state, computational cache and counters
are compared before a checkpoint can be removed. No models are fitted.
}
\details{
With \code{dry_run = FALSE}, the function acquires the mechanism run lock and
per-fit sampler locks, records hashes before deletion, rechecks the retained fit
and every checkpoint, and unlinks only the verified chain RDS files.
Measurement checkpoints, continuation segments, missing, mismatched, ambiguous
and locked fits are preserved. Identity files, timing logs, full fits and all
other outputs remain. Do not run fitting or another cleanup against this output
root concurrently. Existing locks are never overridden.

The manifest reports preserved, verified, removed or partially removed entries.
Receipts are written before deletion and updated afterwards, so an interrupted
cleanup can be audited. Use a fresh audit directory when retrying. This utility
supports the Study II CUT layout and seed convention; it is not a general-purpose
recursive deletion command. Automatic retirement after fit publication uses the
same exact pair verifier, and also conservatively retains continuation segments.
}
\value{A data frame of directory statuses, reasons, matched fits, verified bytes,
removed bytes and number of removed chains.}
\examples{
\dontrun{
prune_study_checkpoints("/path/to/reference",
  audit_dir = "~/checkpoint-audits/dry-run-01")
prune_study_checkpoints("/path/to/reference", dry_run = FALSE,
  audit_dir = "~/checkpoint-audits/cleanup-01")
}
}
