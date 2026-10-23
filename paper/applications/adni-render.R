# Author: Xinye Chen. Fresh-process plotting avoids macOS font initialization after fork.
a<-commandArgs(TRUE);stopifnot(length(a)==5L)
B<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
lib<-Sys.getenv('EIVGP_R_LIBRARY','');if(nzchar(lib)).libPaths(c(lib,.libPaths()))
suppressPackageStartupMessages(library(ggplot2))
p<-readRDS(a[1]);ggsave(a[2],p,width=as.numeric(a[3]),height=as.numeric(a[4]),dpi=as.numeric(a[5]),bg='white')
stopifnot(file.exists(a[2]),file.info(a[2])$size>1000)
