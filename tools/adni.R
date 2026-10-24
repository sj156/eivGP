# ADNI entrypoint; all functions are defined in adni-tools.R.
local({
  files<-as.character(unlist(lapply(sys.frames(),function(frame) {
    c(if(is.character(frame$ofile))frame$ofile,if(is.character(frame$file))frame$file)
  }),use.names=FALSE))
  files<-files[basename(files)=='adni.R']
  path<-if(length(files))tail(files,1L)else sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  sys.source(file.path(dirname(normalizePath(path)),'adni-tools.R'),envir=parent.env(environment()))
})
if(sys.nframe()==0L)adni_method_main('report')
