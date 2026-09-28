## Rscript --vanilla codes/cli/run_competitors.R study1 [plan|run|retry|export] [publication|development]
## Run install_eivgp_dependencies.R explicitly before this module. Never runs MCMC.
args <- commandArgs(trailingOnly=TRUE)
study <- match.arg(if(length(args)) args[1L] else "study1", c("study1","study2"))
action <- match.arg(if(length(args)>1L) args[2L] else "plan", c("plan","run","retry","export"))
mode <- match.arg(if(length(args)>2L) args[3L] else "publication", c("publication","development"))
cli <- grep("^--file=",commandArgs(FALSE),value=TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=","",cli[1L])),"..", ".."),mustWork=TRUE)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
source(file.path(repo,"codes","core/00_diagnostics.R"))
source(file.path(repo,"codes","simulation_helpers.R"))
engine <- mixedgp_simulation_engine(file.path(repo,"codes"))
workers <- suppressWarnings(as.integer(Sys.getenv("EIVGP_WORKERS","1")))
if(length(workers)!=1L || is.na(workers) || workers<1L) stop("EIVGP_WORKERS must be positive.")
data_root <- Sys.getenv("EIVGP_DATA_ROOT",file.path(repo,"reproduction","data","synthetic",study))
cfg <- get(paste0(study,"_simulation_config"))(mode,code_dir=file.path(repo,"codes"),data_root=data_root,workers=workers)
cfg$parallel$workers <- workers
methods <- strsplit(Sys.getenv("EIVGP_COMPETITOR_METHODS","UC-GP,LVGP,EzGP"),",",fixed=TRUE)[[1L]]
methods <- trimws(methods)
if(anyDuplicated(methods) || any(!methods %in% c("UC-GP","LVGP","EzGP"))) stop("Invalid EIVGP_COMPETITOR_METHODS.")
controls <- engine$mixedgp_competitor_protocol(study)
output <- file.path(Sys.getenv("EIVGP_COMPETITOR_REPORT_ROOT",
  file.path(repo,"reproduction","competitor-reports")),study,mode)
overleaf <- Sys.getenv("EIVGP_OVERLEAF_ROOT","")
print(list(study=study,selected_datasets=mode,cells=mixedgp_cell_summary(cfg),
  data_root=data_root,cache=engine$.mixedgp_competitor_cache_root,
  workers=workers,methods=methods,controls=controls,output=output,overleaf=overleaf))
cat("Competitor settings are shared across modes. LVGP's existing time limit is PER ATTEMPT.\n")
if(action=="plan") {
  print(mixedgp_existing_data_status(cfg)); quit(status=0L)
}
source(file.path(repo, "codes", "simulations", "competitor_workflow.R"))
answer <- mixedgp_run_competitor_workflow(cfg, output, workers, methods, action)
if (nzchar(overleaf)) {
  if (!dir.exists(overleaf)) stop("EIVGP_OVERLEAF_ROOT does not exist.")
  target <- file.path(overleaf,"tables","competitors",study,mode)
  dir.create(target,recursive=TRUE,showWarnings=FALSE)
  for (f in list.files(output,full.names=TRUE)) if (!file.copy(f,target,overwrite=TRUE)) stop("Export failed: ",f)
}
cat("Saved competitor-only tables to", output, "\nNo EIV-GP fits were run.\n")
if (any(answer$status$status != "success")) quit(status=1L)
