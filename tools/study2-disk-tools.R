# Study II-specific disk-continuation orchestration; sampler implementation is
# provided by eivGP >=0.4.0.9020. No full-history reads in the sampling executor.
study2_disk_atomic <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile('.commit-',tmpdir=dirname(path));on.exit(unlink(tmp))
  saveRDS(x,tmp,version=3);if(!file.rename(tmp,path))stop('Cannot publish ',path)
}
study2_disk_public <- function(fit,ns=asNamespace('eivGP')) {
  F<-function(n)get(n,ns,inherits=FALSE)
  if(inherits(fit,'eivgp_fit'))return(fit)
  d<-fit$data;cp<-fit$checkpoint;chains<-fit$mcmc$samples_by_chain
  retained<-cp$iteration-cp$burn;nc<-length(cp$states)
  if(!identical(fit$sampler_version,'0.3.1')||!identical(cp$version,2L)||
     !identical(cp$thin,1L)||!identical(cp$control,fit$control)||
     !identical(F('mixedgp_fit_engine')(fit),'multivariate')||
     isTRUE(d$standardize_U)||nc!=4L||length(chains$U)!=nc||
     any(vapply(chains$U,function(x)dim(x)[1L],integer(1))!=retained)||
     any(vapply(cp$states,function(x)is.null(x$rng)||is.null(x$state),logical(1))))
    stop('Unsupported internal Study II checkpoint; no conversion performed')
  if(any(vapply(chains,function(group)length(group)!=nc||any(vapply(group,function(x) {
      if(is.null(x))return(FALSE)
      count<-if(is.null(dim(x)))length(x)else dim(x)[1L]
      count!=retained
    },logical(1))),logical(1))))stop('Retained draw fields have inconsistent lengths')
  if(!is.null(cp$signature)&&!identical(cp$signature,F('mixedgp_v030_checkpoint_hash')(fit)))stop('Modified signed internal checkpoint')
  # Only attach the public continuation metadata that the internal study driver
  # omitted. Data, priors, draws, controls, terminal states and RNG are unchanged.
  before<-list(data=fit$data,priors=fit$priors,control=fit$control,mcmc=fit$mcmc,states=cp$states)
  priors<-fit$priors;priors$dictionary_type<-NULL
  cp$arguments<-list(X_raw=d$X_raw,y_raw=d$y_raw,C_ord=d$C_ord,U_obs=d$U_obs_raw,
    calib_idx=d$calib_idx,d=d$d,m_vec=d$m_vec,ordinal_levels=d$C_level_maps,
    ident=d$ident,kernel=fit$kernel$name,matern_nu=fit$kernel$matern_nu,
    standardize_U=FALSE,store_scores=any(!vapply(chains$S,is.null,logical(1))),
    priors=priors,sampler_control=fit$control,n_iter=cp$iteration,burn=cp$burn,
    thin=1L,n_chains=nc,seed=as.integer(fit$mcmc$chain_stats$seed[1L]-10000L),
    parallel_chains=TRUE,n_cores=4L,verbose=FALSE)
  fit$checkpoint<-cp
  fit$interface<-list(engine='multivariate',parallel=TRUE,n_cores=4L,ordinal_level_maps=d$C_level_maps)
  fit$model_specification<-list(latent_population='N_d(0,I_d)',loading_structure=d$ident,
    ordinal_residual_covariance='I_q',migration='Study II internal checkpoint metadata only')
  class(fit)<-unique(c('eivgp_fit',class(fit)))
  stopifnot(identical(before,list(data=fit$data,priors=fit$priors,control=fit$control,
    mcmc=fit$mcmc,states=fit$checkpoint$states)))
  fit$checkpoint$signature<-F('mixedgp_v030_checkpoint_hash')(fit)
  fit
}
study2_disk_sources <- function(repo,pilot,screen,kind) {
  if(kind=='pilot') {
    manifest<-file.path(pilot,'plan-pilot.csv')
    if(!file.exists(manifest))stop('No pilot manifest at ',manifest,
      '; supply --pilot-root=PATH or EIVGP_STUDY2_PILOT_ROOT')
    keys<-read.csv(manifest,stringsAsFactors=FALSE)
    paths<-list.files(pilot,'^state.rds$',recursive=TRUE,full.names=TRUE)
    if(!length(paths))stop('No saved continuation states')
    jobs<-do.call(rbind,lapply(paths,function(path) {
      id<-basename(dirname(path));rep<-as.integer(sub('rep([0-9]+)-.*','\\1',id))
      cal<-as.integer(sub('.*-cal([0-9]+)$','\\1',id));key<-keys[keys$rep==rep&keys$calibration==cal,]
      if(nrow(key)!=1L)stop('Ambiguous pilot data identity')
      data.frame(job=id,rep=rep,calibration=cal,source=normalizePath(path),data=key$data,stringsAsFactors=FALSE)
    }))
  } else {
    spec<-readRDS(file.path(screen,'experiment.rds'))
    jobs<-data.frame(job=paste0('rep005-cal00-',spec$methods),rep=5L,calibration=0L,
      source=file.path(screen,spec$methods,'fit.rds'),data=spec$data_source,stringsAsFactors=FALSE)
    if(any(!file.exists(jobs$source)))stop('Complete sampler screen fits are required')
  }
  jobs$source_md5<-unname(tools::md5sum(jobs$source))
  jobs$data_md5<-unname(tools::md5sum(jobs$data));jobs
}
study2_disk_main <- function(repo,args=commandArgs(TRUE)) {
  repo<-normalizePath(repo,mustWork=TRUE)
  a<-list(action='plan',source='pilot',cores='12',chunk='500',add='0',round='1',jobs='all',
    output=file.path(repo,'results-sensitivity/study2-disk-continuation-v1'),library='',
    `pilot-root`=Sys.getenv('EIVGP_STUDY2_PILOT_ROOT',file.path(repo,'results-sensitivity/study2-pilot-12cores')),
    `screen-root`=file.path(repo,'results-sensitivity/study2-sampler-screen-4000-v1'),
    `recover-stale-lock`='false')
  if('--help'%in%args) {
    cat(paste('Study II disk-backed continuation',
      'Actions: plan, prepare, migrate, run, diagnose, report. Default: plan; no fitting.',
      '--source=pilot|screen --pilot-root=PATH --screen-root=PATH --output=PATH --library=PATH --jobs=all|job1,job2',
      '--cores=12 --chunk=500 --add=2000 --round=1 --recover-stale-lock=true',
      'prepare pins inputs/runtime; migrate imports serially with n_iter=0.',
      'run requires migrated archives, advances one explicit bounded request, and never loads history.',
      'Repeating the same round resumes its fixed target; a new round opts into another bounded extension.',
      'diagnose is a separate serial full-history analysis command; it can use substantial RAM.',
      'Keep fit.rds and fit.rds.draws together. No original result or frozen pilot is overwritten.',sep='\n'),'\n')
    return(invisible(NULL))
  }
  seen<-character()
  for(x in args) {
    if(!grepl('^--[^=]+=.+$',x))stop('Use --name=value')
    k<-sub('^--([^=]+)=.*$','\\1',x)
    if(!k%in%names(a)||k%in%seen)stop('Unknown or duplicate option: ',k)
    a[[k]]<-sub('^--[^=]+=','',x);seen<-c(seen,k)
  }
  if(!a$action%in%c('plan','prepare','migrate','run','diagnose','report')||!a$source%in%c('pilot','screen'))stop('Invalid action/source')
  for(k in c('cores','chunk','add','round')) {
    if(!grepl('^[0-9]+$',a[[k]]))stop('Invalid integer: ',k)
    a[[k]]<-as.integer(a[[k]]);if(is.na(a[[k]]))stop('Integer overflow')
  }
  if(a$cores<4L||a$chunk<1L||a$round<1L)stop('Need cores>=4, chunk>=1 and round>=1')
  if(!a[['recover-stale-lock']]%in%c('true','false'))stop('Invalid recovery option')
  canonical<-function(p) {
    p<-path.expand(p);if(!grepl('^/',p))p<-file.path(repo,p)
    if(file.exists(p))return(normalizePath(p,mustWork=TRUE))
    file.path(canonical(dirname(p)),basename(p))
  }
  root<-canonical(a$output);pilot<-canonical(a[['pilot-root']]);screen<-canonical(a[['screen-root']])
  overlap<-function(p,q)p==q||startsWith(p,paste0(q,'/'))||startsWith(q,paste0(p,'/'))
  if(any(vapply(c(canonical(file.path(repo,'results')),pilot,screen,file.path(repo,'tools'),file.path(repo,'eivGP')),
    function(p)overlap(root,p),logical(1))))stop('Output must be separate from original results, frozen runs and source')
  specfile<-file.path(root,'plan.rds')
  if(file.exists(specfile)) {
    spec<-readRDS(specfile)
    if(spec$source!=a$source)stop('Prepared source differs; use a separate output root')
    jobs<-spec$jobs
  } else jobs<-study2_disk_sources(repo,pilot,screen,a$source)
  if(a$jobs!='all') {
    ids<-strsplit(a$jobs,',',fixed=TRUE)[[1]]
    if(any(!ids%in%jobs$job)||anyDuplicated(ids))stop('Unknown/duplicate jobs')
    jobs<-jobs[match(ids,jobs$job),,drop=FALSE]
  }
  path_for<-function(job)file.path(root,'jobs',job,'fit.rds')
  if(a$action=='plan') {
    print(jobs[c('job','rep','calibration')],row.names=FALSE)
    cat('Disk buffering:',a$chunk,'draws per chain; four chains/fit; up to',a$cores%/%4L,'concurrent fits.\n')
    cat('Sampling reads only disk indexes, terminal states, and bounded new blocks. Diagnostics are separate.\n')
    cat('Output:',root,'\nNo migration, sampling, or original-result updates.\n');return(invisible(jobs))
  }
  if(a$action!='prepare'&&!file.exists(specfile))stop('Run --action=prepare first')
  lib<-a$library
  if(!nzchar(lib))lib<-if(file.exists(specfile))spec$library else if(dir.exists(file.path(root,'library')))file.path(root,'library')else Sys.getenv('EIVGP_R_LIBRARY',file.path(repo,'.r-library'))
  .libPaths(unique(c(canonical(lib),file.path(repo,'.r-library'),.libPaths())))
  if(!requireNamespace('eivGP',quietly=TRUE)||utils::packageVersion('eivGP')<'0.4.0.9020'||
     !'continue_eivgp_disk'%in%getNamespaceExports('eivGP'))stop('Use an isolated eivGP >=0.4.0.9020 library with --library=PATH')
  F<-function(n)getFromNamespace(n,'eivGP')
  modules<-F('eivgp_module_identity')()
  deps<-c('Rfast','Rcpp','RcppArmadillo','posterior','TruncatedNormal')
  versions<-setNames(vapply(deps,function(n)as.character(utils::packageVersion(n)),character(1)),deps)
  if(file.exists(specfile)&&(!identical(spec$modules,modules)||!identical(spec$dependency_versions,versions)||!identical(spec$R,R.version.string)))stop('Runtime module identity changed; use a new output root')
  if(a$action=='prepare') {
    if(file.exists(specfile))stop('Plan already prepared; use migrate/run/report')
    dir.create(root,recursive=TRUE,showWarnings=FALSE)
    spec<-list(schema=1L,source=a$source,jobs=jobs,modules=modules,dependency_versions=versions,
      library=dirname(find.package('eivGP')),R=R.version.string,created=Sys.time(),
      chunk=a$chunk,diagnostics='Separate, serial, full retained history; never called by sampling')
    study2_disk_atomic(spec,specfile);cat('Prepared immutable input/runtime plan. No sampling.\n');return(invisible(spec))
  }
  index_for<-function(job)F('mixedgp_disk_index')(path_for(job))
  if(a$action=='report') {
    rows<-lapply(jobs$job,function(job) {
      p<-path_for(job)
      if(!file.exists(file.path(paste0(p,'.draws'),'index.rds')))return(data.frame(job=job,iteration=NA,draw_blocks=0L,index_kib=NA,status='not_migrated'))
      z<-index_for(job);data.frame(job=job,iteration=z$checkpoint$iteration,draw_blocks=nrow(z$disk_history$blocks),
        index_kib=as.numeric(object.size(z))/1024,status='disk_checkpoint')
    });print(do.call(rbind,rows),row.names=FALSE);return(invisible(rows))
  }
  ps<-suppressWarnings(system2('ps',c('-axo','pid=,command='),stdout=TRUE,stderr=TRUE))
  if(!is.null(attr(ps,'status')))stop('Cannot inspect process ownership safely')
  if(any(grepl('study2-targeted.R',ps,fixed=TRUE)&grepl(pilot,ps,fixed=TRUE))||
     any(grepl('study2-sampler-pilot.R',ps,fixed=TRUE)&grepl('--action=run',ps,fixed=TRUE)))
    stop('Another Study II sampler is active; avoid overlapping CPU/memory budgets')
  pids<-as.integer(sub('^ *([0-9]+).*','\\1',ps))
  if(any(pids!=Sys.getpid() & grepl('/exec/R ',ps,fixed=TRUE) &
    grepl('study2-continue(-eivGP)?[.][Rr]',ps) & grepl('--action=run',ps,fixed=TRUE)))
    stop('Another disk sampling phase is active; migration and full-history analysis must wait')
  lock<-file.path(root,'LOCK')
  if(dir.exists(lock)) {
    owner<-file.path(lock,'owner.rds')
    if(!file.exists(owner))stop('Unknown lock ownership')
    own<-readRDS(owner)
    if(own$host!=Sys.info()[['nodename']]||own$pid%in%pids)stop('Disk continuation is already owned by an active or remote process')
    if(a[['recover-stale-lock']]!='true')stop('Stale app lock; inspect processes then use --recover-stale-lock=true')
    unlink(lock,recursive=TRUE)
  }
  if(!dir.create(lock,showWarnings=FALSE))stop('Another writer acquired the archive lock')
  study2_disk_atomic(list(pid=Sys.getpid(),host=Sys.info()[['nodename']]),file.path(lock,'owner.rds'))
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  if(!identical(unname(tools::md5sum(jobs$source)),jobs$source_md5)||
     !identical(unname(tools::md5sum(jobs$data)),jobs$data_md5))stop('Source input hash changed')
  if(a$action=='migrate') {
    # Deliberately serial: import can read one legacy fit into RAM, but never
    # holds multiple histories before forking continuation workers.
    for(i in seq_len(nrow(jobs))) {
      job<-jobs[i,];path<-path_for(job$job)
      if(file.exists(file.path(paste0(path,'.draws'),'index.rds'))) {F('mixedgp_disk_index')(path,verify=TRUE);next}
      if(!file.exists(path)) {
        source<-F('mixedgp_read_artifact')(job$source)
        fit<-source[['fit',exact=TRUE]];if(is.null(fit))fit<-source
        frozen<-readRDS(job$data);tr<-frozen$data$train
        eq<-function(x,y)isTRUE(all.equal(unname(x),unname(y),check.attributes=FALSE,tolerance=1e-12))
        expected<-matrix(NA_real_,nrow(tr$U),ncol(tr$U));ids<-frozen$calibration_sets[[as.character(job$calibration)]]
        expected[ids,]<-tr$U[ids,,drop=FALSE]
        if(!eq(fit$data$X_raw,tr$X)||!eq(fit$data$y_raw,tr$y)||!eq(fit$data$C_ord,tr$C)||!eq(fit$data$U_obs,expected))stop('Frozen data mismatch: ',job$job)
        fit<-study2_disk_public(fit)
        # Atomic serialization only: do not call save_fit on legacy objects,
        # which could retire source-referenced redundant restart files.
        study2_disk_atomic(F('mixedgp_pack_artifact')(fit),path)
        rm(source,fit,frozen,tr);gc(FALSE)
      }
      eivGP::continue_eivgp_disk(path,0L,chunk_iter=a$chunk,parallel=FALSE,n_cores=1L)
      gc(FALSE);cat('Imported without sampling:',job$job,'\n')
    }
    return(invisible(NULL))
  }
  if(any(!file.exists(vapply(jobs$job,function(j)file.path(paste0(path_for(j),'.draws'),'index.rds'),character(1)))))
    stop('Run migrate first; sampling never imports full histories')
  if(a$action=='diagnose') {
    for(i in seq_len(nrow(jobs))) {
      job<-jobs[i,];path<-path_for(job$job);index<-index_for(job$job)
      message('Full-history analysis, serial and separate from MCMC: ',job$job)
      fit<-eivGP::read_fit(path,fields=setdiff(index$disk_history$fields,'S'))
      test<-readRDS(job$data)$data$test;ids<-F('mixedgp_study2_diagnostic_rows')(test$X,test$C,5L)
      series<-list(scientific=F('mixedgp_continuation_target_series')(fit,test$X[ids,,drop=FALSE],test$C[ids,,drop=FALSE],
        if(job$calibration>0)test$U[ids,,drop=FALSE]else NULL,64L,1000000L+10000L*job$rep+2900L,TRUE,4L),
        invariant=F('mixedgp_study2_invariant_series')(fit),raw=F('mixedgp_study2_raw_series')(fit))
      tabs<-lapply(series,F('mixedgp_summarize_diagnostic_series'))
      folder<-file.path(dirname(path),paste0('analysis-',index$checkpoint$iteration));dir.create(folder,showWarnings=FALSE)
      for(k in names(tabs))write.csv(tabs[[k]],file.path(folder,paste0(k,'.csv')),row.names=FALSE)
      required<-c('scientific','invariant',if(job$calibration>0)'raw')
      pass<-all(vapply(required,function(k)F('mixedgp_diagnostic_table_pass')(tabs[[k]],1.01,400,400,expected_parameters=names(series[[k]])),logical(1)))
      table<-do.call(rbind,tabs[required]);pass<-pass&&all(is.finite(table$mcse_sd_ratio)&table$mcse_sd_ratio<=.05)
      study2_disk_atomic(list(iteration=index$checkpoint$iteration,full_retained_history=TRUE,pass=pass),file.path(folder,'assessment.rds'))
      rm(fit,series,tabs,table,index);gc(FALSE)
    };return(invisible(NULL))
  }
  if(a$add<1L)stop('Specify a positive --add for a bounded sampling request; no automatic extensions')
  available<-parallel::detectCores(logical=TRUE)
  if(is.na(available)||available<4L)stop('Hardware detection unavailable; run in native terminal')
  workers<-min(nrow(jobs),min(a$cores,available)%/%4L)
  Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  run_one<-function(job)tryCatch({
    path<-path_for(job);index<-index_for(job)
    request<-file.path(dirname(path),sprintf('request-%03d.rds',a$round))
    if(file.exists(request)) {
      target<-readRDS(request)
      if(target$add!=a$add||target$chunk!=a$chunk)stop('Request changed; use the original add/chunk when resuming')
    } else {
      if(a$round>1L) {
        previous<-file.path(dirname(path),sprintf('request-%03d.rds',a$round-1L))
        if(!file.exists(previous)||readRDS(previous)$to!=index$checkpoint$iteration)stop('Previous round must be complete')
      }
      target<-list(from=index$checkpoint$iteration,to=index$checkpoint$iteration+a$add,add=a$add,chunk=a$chunk)
      study2_disk_atomic(target,request)
    }
    if(index$checkpoint$iteration>target$to)stop('Archive already advanced beyond this request')
    progress<-function(z) {
      stopifnot(is.null(z$mcmc$samples_by_chain))
      row<-data.frame(job=job,iteration=z$checkpoint$iteration,target=target$to,chunk=a$chunk,
        index_kib=as.numeric(object.size(z))/1024,updated_utc=format(Sys.time(),tz='UTC',usetz=TRUE))
      write.csv(row,file.path(dirname(path),'progress.csv'),row.names=FALSE)
    }
    eivGP::continue_eivgp_disk(path,target$to-index$checkpoint$iteration,chunk_iter=a$chunk,
      parallel=TRUE,n_cores=4L,verbose=TRUE,progress=progress)
    list(ok=TRUE)
  },error=function(e){writeLines(conditionMessage(e),file.path(dirname(path_for(job)),'error.txt'));list(ok=FALSE,error=conditionMessage(e))})
  results<-parallel::mclapply(jobs$job,run_one,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  if(any(!vapply(results,function(x)is.list(x)&&isTRUE(x$ok),logical(1))))stop('Some disk requests failed; inspect job error.txt files')
  cat('Bounded sampling finished. No history materialization or diagnostic assessment was performed.\n')
  invisible(results)
}
