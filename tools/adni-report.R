# ADNI command-line entrypoint.
if(sys.nframe()==0L) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  source(file.path(dirname(normalizePath(script)),'adni-tools.R'))
  adni_method_main('report')
}
