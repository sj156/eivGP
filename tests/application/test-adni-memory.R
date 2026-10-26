repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/adni.R'))
list2env(as.list(adni_load_api(repo)),envir=environment())

testthat::test_that('ADNI ESS settings are validated',{
  expect_null(adni_settings()$ess_min)
  expect_equal(adni_settings(ess_min='100')$ess_min,100)
  expect_equal(adni_cli_params('--ess-min=100')$ess_min,'100')
  expect_error(adni_settings(ess_min='nonsense'),'ess_min')
  m<-list(stages=c(100L,200L),burn=20L,folds=1L,cells='O50',ess_min=100)
  expect_null(adni_screen_settings(m,100L,1L)$ess_min)
  expect_equal(adni_screen_settings(m,200L,1L)$ess_min,100)
})

testthat::test_that('ADNI saved fit retains scientific diagnostics',{
  set.seed(19);X<-matrix(rnorm(24),12,2);U<-matrix(rnorm(24),12,2)
  C<-matrix(sample(1:4,48,TRUE),12,4);y<-X[,1]+U[,1]+rnorm(12)
  fit<-eivGP::fit_eivgp(X,y,C,U_obs=U,engine='multivariate',latent_dim=2,
    m_vec=rep(4L,4),ident='lower_triangular',standardize_U=FALSE,
    n_iter=22,burn=10,n_chains=4,parallel=FALSE,seed=7)
  # Actual ADNI task routes a saved fit through final-stage science checks.
  out<-tempfile();dir.create(out);on.exit(unlink(out,recursive=TRUE))
  task<-data.frame(fold=1L,cell='O50')
  path<-adni_plan_fit_path(out,task);adni_plan_save_fit(fit,path)
  te<-data.frame(X1=X[,1],X2=X[,2],C,U1=U[,1],U2=U[,2]);names(te)[3:6]<-paste0('C',1:4)
  z<-list(test=te,x=c('X1','X2'),c=paste0('C',1:4),u=c('U1','U2'),protocol=list(seeds=list(prediction=11L)))
  s<-list(iter=22L,burn=10L,chains=4L,cores=1L,final_stage=TRUE,ess_min=.01,samplers=list(O50=fit$control$sampler))
  suppressWarnings(adni_plan_fit_task(task,list(cells=list('1/O50'=z)),out,s))
  saved<-adni_plan_read_fit(path)
  expect_true(saved$scientific_ess$success)
  expect_identical(saved$mcmc,fit$mcmc)
  expect_true(file.exists(file.path(dirname(path),'scientific-trace-targets.rds')))
})

testthat::test_that('legacy ADNI category encoding is preserved',{
expect_identical(adni_make_C(data.frame(c_ab42_ab40_3=c(1,3),c_gfap_3=c(2,1)))$c_gfap_3, ordered(c(2,1),levels=1:3))
})

testthat::test_that('one core budget determines both phase allocations',{
 for(cores in 1:16)for(chains in c(1L,4L,8L))for(tasks in c(0L,1L,15L)) {
  expected<-as.integer(min(tasks,max(1L,floor(cores/chains))))
  for(phase in c('fit','evaluate'))testthat::expect_identical(
    adni_resident_workers(cores,phase,tasks,chains),expected)
 }
 testthat::expect_identical(adni_resident_workers(6,'fit',15,4),1L)
 for(v in c('0','-1','1.5','bad'))testthat::expect_error(adni_resident_workers(v,'fit',15,4),'integer')
 testthat::expect_identical(adni_cli_params('--cores=12'),list(cores='12'))
 for(flag in c('--fit-jobs=1','--evaluation-jobs=1'))testthat::expect_error(adni_cli_params(flag),'Unknown')
 testthat::expect_error(adni_method_run(repo,list(fit_jobs=1),'eivGP'),'only --cores')
})

testthat::test_that('continuation reads a checkpoint once and preserves validation',{
 e<-adni_load_api(repo);path<-tempfile();file.create(path);on.exit(unlink(path))
 reads<-0L;e$adni_plan_fit_path<-function(...)path
 e$adni_plan_read_fit<-function(...) {reads<<-reads+1L;list(checkpoint=list(iteration=10L,burn=2L,states=list(1,2,3,4)),control=list(sampler='ESS'))}
 e$adni_plan_diagnostics<-function(...)list(parameters=data.frame(),summary=data.frame())
 e$mixedgp_atomic_csv<-function(...)NULL
 task<-list(fold=1,cell='O50');s<-list(iter=10L,burn=2L,chains=4L,samplers=list(O50='ESS'))
 e$adni_plan_fit_task(task,list(cells=list('1/O50'=list())),tempdir(),s)
 testthat::expect_identical(reads,1L)
 s$burn<-3L;testthat::expect_error(e$adni_plan_fit_task(task,list(cells=list()),tempdir(),s),'sampler/burn')
 s$burn<-2L;s$iter<-9L;testthat::expect_error(e$adni_plan_fit_task(task,list(cells=list()),tempdir(),s),'earlier stage')
})
