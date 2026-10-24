# Ocean entrypoint; all functions are defined in ocean-tools.R.
local({
  files<-as.character(unlist(lapply(sys.frames(),function(frame) {
    c(if(is.character(frame$ofile))frame$ofile,if(is.character(frame$file))frame$file)
  }),use.names=FALSE))
  files<-files[basename(files)=='ocean.R']
  path<-if(length(files))tail(files,1L)else sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  sys.source(file.path(dirname(normalizePath(path)),'ocean-tools.R'),envir=parent.env(environment()))
})
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(identical(args,"--help")) {
    cat("Rscript tools/ocean.R --action=plan|fit|analyze|continue|report --profile=paper|pilot|quick\n",
      "  --dataset=ocean200|ocean500|ocean1000 --cores=16 --output=/path/to/results\n",
      "  Default output: results/<dataset>/ (paper), results/<dataset>_<profile>/ otherwise.\n",
      "  Optional: --data-root=... --folds=1:5 --cells=O0,O20,O50 --modules=...\n",
      "  --run-id=... --diagnostics-dir=...\n",
      "  Sampling is fixed by profile: paper 4 chains/50000 iterations/5000 burn-in; pilot 4/80/30; quick 1/80/30.\n",
      "  Incompatible initial-fit results are preserved; a new folder is selected and printed.\n",
      "  --continuation-iter=20000 (additional iterations per existing chain; pilot default: 80)\n",
      "  --task-ids=... --import-dirs=... --cc-output-dir=... --representative-fold=1\n",
      "No arguments: print the plan. --cores=N alone: fit, analyze, export paper assets and render ocean500.\n",
      "fit (also run) fits, analyzes with the shared core budget and writes HTML; continue extends, analyzes and updates HTML.\n",
      "cores is the shared budget across fitting, prediction, diagnostics and plotting. analyze/report rebuild HTML without fitting. See README_ocean.md.\n",sep="")
  } else {
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    root<-normalizePath(file.path(dirname(normalizePath(script)),".."))
    Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
    lib<-Sys.getenv("EIVGP_R_LIBRARY",file.path(root,".r-library"))
    if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
    ocean_run_cli(root,ocean_cli_params(args))
  }
}
