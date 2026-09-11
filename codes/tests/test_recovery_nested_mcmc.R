source('codes/simulation_helpers.R');source('codes/publication_recovery.R')
results<-list()
mixedgp_recovery_dispatch(as.list(1:2),function(i){
  set.seed(i)
  sys.source('codes/tests/test_publication_recovery.R',envir=new.env(parent=globalenv()))
  list(ok=TRUE,pid=Sys.getpid())
},function(z,i){results[[i]]<<-z},workers=2L)
stopifnot(length(results)==2L,all(vapply(results,function(z)isTRUE(z$ok),logical(1))),
 length(unique(vapply(results,`[[`,integer(1),'pid')))==2L)
message('Two concurrent real tiny MCMC dataset drivers, each with two forked chains, passed.')
