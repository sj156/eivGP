# Ocean preparation entrypoint; functions are in ocean-tools.R.
if(sys.nframe()==0L) {
  entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
  root<-normalizePath(file.path(dirname(entry),'..'))
  source(file.path(root,'tools/ocean-tools.R'))
  args<-commandArgs(TRUE)
  ocean_prepare_main(args,repository_dir=root)
}
