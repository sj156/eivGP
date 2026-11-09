#!/usr/bin/env Rscript
# EIV-GP continuation entrypoint. Default: plan. Study II continuation uses bounded disk history, not the old
# frozen in-memory runner. Existing fit archives remain read-only inputs.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1L])
if(sys.nframe()==0L) {
  source(file.path(dirname(normalizePath(script)),'study2-disk-tools.R'))
  study2_disk_main(dirname(dirname(normalizePath(script))))
}
