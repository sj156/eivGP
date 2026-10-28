# Author: Xinye Chen. OAI command-line entrypoint.
if(sys.nframe()==0L) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  source(file.path(dirname(normalizePath(script)),'oai-tools.R'))
  oai_method_main('report')
}
