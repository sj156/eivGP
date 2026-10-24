\name{paper-workflows}
\alias{study_config}
\alias{run_numerical_study}
\alias{read_study_results}
\alias{summarize_comparison}
\alias{summarize_diagnostics}
\alias{plot_comparison}
\alias{report_table}
\alias{read_fit}
\alias{save_fit}
\title{Run analyses and report saved results}
\description{Package-owned numerical and application workflows, reusable score summaries, plots, tables, and lossless fit storage. Paper scripts select scientific settings and external data/output locations.}
\usage{
study_config(study = c("study1", "study2"),
  profile = c("quick", "development", "paper"), output_root, data_root, cores = 4L)
run_numerical_study(config, action = c("plan", "data", "fit", "run", "report"),
  competitors = c("UC-GP", "LVGP", "EzGP"), run_dir = NULL)
read_study_results(path, competitor_dir = NULL)
summarize_comparison(results, metrics = c("RMSE", "CRPS", "Coverage95", "Width95"))
summarize_diagnostics(diagnostics)
plot_comparison(summary, metric = "RMSE")
report_table(x, caption = NULL, digits = 3L)
read_fit(path)
save_fit(object, path)
}
\arguments{
\item{study}{Numerical design identifier.}
\item{profile}{Scientific budget profile. Quick runs are workflow examples.}
\item{output_root, data_root, path, run_dir, competitor_dir}{Explicit external directories or artifact paths.}
\item{cores}{Total CPU budget.}
\item{config}{Resolved settings returned by \code{study_config}.}
\item{action}{Planning, execution, or reporting. Report mode does not fit models.}
\item{competitors}{Competitor package adapters to run; an empty character vector skips them.}
\item{results}{A data frame or saved-result list from \code{read_study_results}.}
\item{metrics, metric}{Score names.}
\item{diagnostics}{Saved per-fit diagnostics.}
\item{summary}{Aggregated comparison table.}
\item{x}{Table data.}
\item{caption, digits}{Table presentation settings.}
\item{object}{Fit or container to save losslessly.}
}
\value{Configuration or result lists, summary data frames, a ggplot object, a formatted table, or a restored fit, depending on the function. Saving returns the output path invisibly.}
\details{Existing result archives are never migrated automatically. Compact saves retain chain-resolved histories and terminal states, rebuilding identical pooled views on read. Replication MCSEs do not quantify within-fit MCMC error. Fitting uses existing package implementations and competitor packages; no automatic chain extension or diagnostic filtering is introduced.}
