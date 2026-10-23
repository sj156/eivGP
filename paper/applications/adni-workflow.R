# Author: Xinye Chen. Shared CLI orchestration; implementation is in eivGP.
adni_method_main <- function(group,args=commandArgs(TRUE)) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  repo<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
  loader<-new.env(parent=globalenv())
  sys.source(file.path(repo,'paper/applications/adni.R'),loader)
  api<-loader$adni_load_api(repo)
  if(identical(args,'--help')) {
    cat('Rscript paper/applications/adni-',group,'.R --profile=paper --cores=4\n',
      '  --profile=pilot|paper|quick; --dataset=NAME; --data-root=PATH; --output=PATH\n',
      '  --config=FILE; --iter=N --burn=N --chains=N --stages=N,N\n',
      '  --folds=1:5 --cells=O0,O20,O50; --plan=true; --reuse-only=true; --report-source=FILE\n',
      '  Run EIV-GP and competitors sequentially; report reads saved results.\n',sep='')
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  api$adni_method_run(repo,api$adni_cli_params(args),group)
}
