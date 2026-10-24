# ADNI preparation entrypoint; functions are in adni-tools.R.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(sys.nframe()==0L) {
  source(file.path(dirname(normalizePath(script)),'adni-tools.R'))
  args<-commandArgs(TRUE)
  if(length(args)!=3L)stop('Usage: Rscript tools/adni-prepare-folds.R PRIVATE_COHORT.csv PRIVATE_METADATA.json OUTPUT_INPUT_DIR')
  adni_prepare_folds(args[1],args[2],args[3])
}
