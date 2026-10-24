# Ocean preparation entrypoint; functions are in ocean-tools.R.
if(sys.nframe()==0L) {
  entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
  root<-normalizePath(file.path(dirname(entry),'..'))
  source(file.path(root,'tools/ocean-tools.R'))
  args<-commandArgs(TRUE)
  prepare_ocean_folds(if(length(args))args[1]else'ocean1000',if(length(args)>1L)args[2]else NULL,repository_dir=root)
}
