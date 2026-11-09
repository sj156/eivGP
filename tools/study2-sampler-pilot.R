#!/usr/bin/env Rscript
# Separate rep005/cal00 ESS/HMC/HYBRID experiment. Default: plan, never fitting.
# prepare pins source and installs an isolated package; run explicitly fits.
# Existing primary results and continuation-pilot artifacts are read-only inputs.
study2_sampler_pilot <- function(repo, args=commandArgs(TRUE)) {
  repo<-normalizePath(repo,mustWork=TRUE)
  a<-list(action='plan',cores='12',iterations='4000',warmup='2000',
    output=file.path(repo,'results-sensitivity/study2-sampler-screen-4000-v1'),
    `existing-pilot`=Sys.getenv('EIVGP_STUDY2_PILOT_ROOT',file.path(repo,'results-sensitivity/study2-pilot-12cores')),
    `recover-stale-locks`='false')
  if('--help'%in%args) {
    cat('Actions: plan, prepare, validate, run, report. Options: --output=PATH --cores=12 --iterations=4000 --warmup=2000 --existing-pilot=PATH --recover-stale-locks=true\n')
    cat('prepare requires a fresh output root and installs only its own library. run refuses to overlap the current continuation pilot.\n')
    return(invisible(NULL))
  }
  seen<-character()
  for(x in args) {
    if(!grepl('^--[^=]+=.+$',x))stop('Use --name=value')
    k<-sub('^--([^=]+)=.*$','\\1',x)
    if(!k%in%names(a)||k%in%seen)stop('Unknown or duplicate option: ',k)
    a[[k]]<-sub('^--[^=]+=','',x);seen<-c(seen,k)
  }
  if(!a$action%in%c('plan','prepare','validate','run','report'))stop('Invalid action')
  for(k in c('cores','iterations','warmup')) {
    if(!grepl('^[0-9]+$',a[[k]]))stop(k,' must be an integer')
    a[[k]]<-as.integer(a[[k]]);if(is.na(a[[k]]))stop('Integer overflow')
  }
  if(a$cores<4L||a$warmup<4L||a$iterations<a$warmup+4L)stop('Use cores>=4, warmup>=4, and at least four retained draws')
  if(!a[['recover-stale-locks']]%in%c('true','false'))stop('Invalid recover-stale-locks')
  canonical<-function(p) {
    p<-path.expand(p);if(!grepl('^/',p))p<-file.path(repo,p)
    if(file.exists(p))return(normalizePath(p,mustWork=TRUE))
    file.path(canonical(dirname(p)),basename(p))
  }
  root<-canonical(a$output);existing<-canonical(a[['existing-pilot']])
  overlap<-function(p,q)p==q||startsWith(p,paste0(q,'/'))||startsWith(q,paste0(p,'/'))
  if(overlap(root,canonical(file.path(repo,'results')))||overlap(root,existing)||
     overlap(root,file.path(repo,'eivGP'))||overlap(root,file.path(repo,'tools')))
    stop('Output must be separate from existing results, pilot and source trees')
  atomic<-function(x,p) {
    dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
    tmp<-tempfile('.publish-',tmpdir=dirname(p));on.exit(unlink(tmp))
    saveRDS(x,tmp,version=3);if(!file.rename(tmp,p))stop('Atomic publication failed: ',p)
  }
  .libPaths(unique(c(file.path(repo,'.r-library'),.libPaths())))
  planfile<-file.path(root,'experiment.rds')
  # Reporting reads summaries only; it neither loads fit histories nor fits.
  report<-function() {
    rows<-lapply(c('ESS','HMC','HYBRID'),function(method) {
      p<-file.path(root,method,'assessment.rds')
      if(file.exists(p))readRDS(p)$summary else {
        err<-file.path(root,method,'error.txt')
        data.frame(method=method,status=if(file.exists(err))'failed'else'pending')
      }
    })
    columns<-unique(unlist(lapply(rows,names)))
    rows<-lapply(rows,function(r){for(k in setdiff(columns,names(r)))r[[k]]<-NA;r[columns]})
    tab<-do.call(rbind,rows);print(tab,row.names=FALSE)
    if(dir.exists(root)) {
      write.csv(tab,file.path(root,'sampler-comparison.csv'),row.names=FALSE)
      comparisons<-list()
      baseline<-file.path(root,'ESS/scientific-diagnostics.csv')
      if(file.exists(baseline)&&file.exists(file.path(root,'ESS/assessment.rds')))for(method in c('HMC','HYBRID')) {
        target<-file.path(root,method,'scientific-diagnostics.csv')
        if(!file.exists(target)||!file.exists(file.path(root,method,'assessment.rds')))next
        ess<-read.csv(baseline);other<-read.csv(target)
        both<-merge(ess[,c('parameter','mean','mcse')],other[,c('parameter','mean','mcse')],by='parameter',suffixes=c('_ESS','_alternative'))
        both$method<-method;both$difference<-both$mean_alternative-both$mean_ESS
        both$combined_mcse<-sqrt(both$mcse_ESS^2+both$mcse_alternative^2)
        both$difference_over_mcse<-both$difference/both$combined_mcse
        both$both_methods_passed<-all(tab$status[tab$method%in%c('ESS',method)]=='passed')
        comparisons[[method]]<-both
      }
      if(length(comparisons))write.csv(do.call(rbind,comparisons),file.path(root,'scientific-target-agreement.csv'),row.names=FALSE)
    }
    cat('\nFresh, matched-budget fits only. Missing values are not passes.\n')
    if('min_bulk_per_second'%in%names(tab)) {
      usable<-which(tab$status=='passed'&is.finite(tab$min_bulk_per_second))
      if(length(usable))cat('Highest minimum bulk ESS/s among passing fits in this one pilot: ',tab$method[usable[which.max(tab$min_bulk_per_second[usable])]],'\n',sep='')
      else cat('No passing sampler with comparable timing yet; no winner selected.\n')
    }
    history<-file.path(existing,'rep005-cal00/history.csv')
    if(file.exists(history)) {
      cat('\nHistorical ESS continuation (different budget/source; context only):\n')
      print(read.csv(history)[,c('iteration','max_rhat','min_bulk','min_tail','status')],row.names=FALSE)
    }
    invisible(tab)
  }
  if(a$action=='report')return(report())
  if(a$action%in%c('run','validate')&&!file.exists(planfile))stop('Run prepare first, in a fresh output root')
  if(a$action%in%c('plan','prepare')) {
    manifest<-file.path(existing,'plan-pilot.csv')
    if(!file.exists(manifest))stop('No pilot manifest at ',manifest,
      '; supply --existing-pilot=PATH or EIVGP_STUDY2_PILOT_ROOT')
    selection<-read.csv(manifest,stringsAsFactors=FALSE)
    selection<-selection[selection$rep==5L&selection$calibration==0L,]
    if(nrow(selection)!=1L)stop('Need exactly one archived rep005/cal00 reference')
    archive<-readRDS(selection$fit);ref<-archive
    # Exact indexing avoids partial matching of compact-fit mcmc_names.
    while(is.null(ref[['mcmc',exact=TRUE]])&&!is.null(ref[['fit',exact=TRUE]]))ref<-ref[['fit']]
    if(!is.list(ref[['mcmc',exact=TRUE]]))stop('Unsupported reference archive')
    frozen<-readRDS(selection$data);train<-frozen$data$train
    equal<-function(x,y)isTRUE(all.equal(unname(x),unname(y),check.attributes=FALSE,tolerance=1e-12))
    if(!equal(ref$data$X_raw,train$X)||!equal(ref$data$y_raw,train$y)||!equal(ref$data$C_ord,train$C)||any(is.finite(ref$data$U_obs)))
      stop('Archived data/calibration mismatch')
    controls<-lapply(c('ESS','HMC','HYBRID'),function(method) {
      c<-ref$control;c$sampler<-method
      if(method=='HYBRID')c$hybrid_every<-20L
      c
    });names(controls)<-c('ESS','HMC','HYBRID')
    p<-list(schema=1L,rep=5L,calibration=0L,methods=names(controls),
      iterations=a$iterations,warmup=a$warmup,n_chains=4L,
      seed=as.integer(ref$mcmc$chain_stats$seed[1L]-10000L),
      fit_source=normalizePath(selection$fit),data_source=normalizePath(selection$data),
      source_md5=unname(tools::md5sum(c(selection$fit,selection$data))),
      normalized_priors=ref$priors,controls=controls,reference_initial=ref$mcmc$chain_initial,
      store_scores=any(!vapply(ref$mcmc$samples_by_chain$S,is.null,logical(1))),
      model_data=ref$data,kernel=ref$kernel,diagnostic_seed=1052900L,
      n_latent=64L,panel_points=5L,original_continuation_root=existing,
      comparison='Fresh matched-budget ESS/HMC/HYBRID; historical ESS is context only',
      acceptance=list(rhat=1.01,bulk=400,tail=400,mcse_sd=.05),
      timing='ESS/s only for uninterrupted fits; resumed timing is marked unavailable',
      library=file.path(root,'library'),snapshot=file.path(root,'source','eivGP'))
    cat('Rep005/cal00; fixed archived data, priors, model, and four initial chain states.\n')
    print(data.frame(method=p$methods,chains=4L,iterations=p$iterations,warmup=p$warmup,
      hybrid_ESS_sweep_every=c(NA,NA,20L)),row.names=FALSE)
    cat('Seed:',p$seed,'; maximum concurrent fits:',min(3L,a$cores%/%4L),'\nOutput:',root,'\n')
    if(a$action=='plan')return(invisible(p))
    if(dir.exists(root)&&length(list.files(root,all.files=TRUE,no..=TRUE)))stop('prepare requires a fresh output root')
    dir.create(file.path(root,'source'),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(file.path(repo,'eivGP'),file.path(root,'source'),recursive=TRUE))stop('Cannot snapshot package')
    modules<-c(file.path(p$snapshot,'R',c('helpers.R','mcmc.R','hmc.R','RcppExports.R','plotting.R','disk.R')),
      list.files(file.path(p$snapshot,'src'),'([.]cpp$|^Makevars)',full.names=TRUE))
    write.csv(data.frame(module=basename(modules),md5=unname(tools::md5sum(modules))),file.path(p$snapshot,'inst/module-manifest.csv'),row.names=FALSE)
    # Include all R/native source, not merely the module manifest, in identity.
    files<-c(list.files(file.path(p$snapshot,'R'),full.names=TRUE),
      list.files(file.path(p$snapshot,'src'),'([.]cpp$|^Makevars)',full.names=TRUE),
      file.path(p$snapshot,c('DESCRIPTION','NAMESPACE','inst/module-manifest.csv')))
    originals<-file.path(repo,'eivGP',substring(files,nchar(p$snapshot)+2L))
    scientific<-!grepl('inst/module-manifest.csv$',files)
    if(!identical(unname(tools::md5sum(files[scientific])),unname(tools::md5sum(originals[scientific]))))stop('Source changed during snapshot; prepare a fresh root')
    p$snapshot_hashes<-setNames(unname(tools::md5sum(files)),substring(files,nchar(p$snapshot)+2L))
    file.copy(file.path(repo,'renv.lock'),file.path(root,'source','renv.lock'))
    dir.create(p$library,recursive=TRUE,showWarnings=FALSE)
    atomic(p,planfile)
    Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
    log<-file.path(root,'prepare-install.log')
    code<-system2(file.path(R.home('bin'),'R'),c('CMD','INSTALL','--preclean',shQuote(paste0('--library=',p$library)),shQuote(p$snapshot)),stdout=log,stderr=log)
    if(code!=0L)stop('Isolated installation failed; see ',log)
    cat('Isolated package prepared. No scientific fitting was launched.\n')
    return(invisible(p))
  }
  p<-readRDS(planfile)
  if(p$iterations!=a$iterations||p$warmup!=a$warmup)stop('Requested budget differs from prepared plan')
  if(a$action=='run') {
  processes<-suppressWarnings(system2('ps',c('-axo','pid=,command='),stdout=TRUE,stderr=TRUE))
  if(!is.null(attr(processes,'status')))stop('Cannot inspect running processes; no fitting launched')
  if(any(grepl('study2-targeted.R',processes,fixed=TRUE)&grepl(existing,processes,fixed=TRUE)))
    stop('Existing 12-core continuation pilot is active. Run this comparison after it stops to avoid oversubscription.')
  available<-parallel::detectCores(logical=TRUE)
  if(is.na(available)||available<4L)stop('Hardware detection unavailable; use the native terminal')
  cores<-min(a$cores,available);workers<-min(3L,cores%/%4L)
  lock<-file.path(root,'RUNNING')
  recovery<-a[['recover-stale-locks']]=='true'
  pids<-as.integer(sub('^ *([0-9]+).*','\\1',processes))
  if(dir.exists(lock)) {
    owner<-file.path(lock,'owner');pid<-if(file.exists(owner))as.integer(readLines(owner)[1])else NA_integer_
    if(!is.na(pid)&&pid%in%pids)stop('Sampler comparison already running')
    if(!recovery)stop('Stale comparison lock; inspect processes, then use --recover-stale-locks=true')
    unlink(lock,recursive=TRUE)
  }
  if(!dir.create(lock,showWarnings=FALSE))stop('Cannot acquire comparison lock')
  writeLines(as.character(Sys.getpid()),file.path(lock,'owner'));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  if(recovery)for(path in list.dirs(root,recursive=TRUE,full.names=TRUE)) {
    if(basename(path)!='LOCK')next
    owner<-file.path(path,'owner.csv')
    if(!file.exists(owner))stop('Unknown checkpoint lock: ',path)
    owners<-read.csv(owner)
    if(any(owners$pid%in%pids)||any(owners$host!=Sys.info()[['nodename']]))stop('Checkpoint owner may be active: ',path)
    unlink(path,recursive=TRUE)
  }
  }
  if(!identical(unname(tools::md5sum(file.path(p$snapshot,names(p$snapshot_hashes)))),unname(p$snapshot_hashes))||
     !identical(unname(tools::md5sum(c(p$fit_source,p$data_source))),p$source_md5))stop('Frozen source/input identity changed')
  .libPaths(unique(c(p$library,.libPaths())))
  if(!requireNamespace('jsonlite',quietly=TRUE))stop('jsonlite is required to verify the dependency lock')
  dependency_lock<-jsonlite::read_json(file.path(root,'source/renv.lock'))
  expected<-vapply(dependency_lock$Packages,`[[`,character(1),'Version')
  expected[]<-as.character(package_version(expected))
  actual<-vapply(names(expected),function(n)if(requireNamespace(n,quietly=TRUE))as.character(utils::packageVersion(n))else 'missing',character(1))
  if(any(actual!=expected))stop('Pinned dependencies unavailable or different: ',paste(names(expected)[actual!=expected],collapse=', '))
  F<-function(n)getFromNamespace(n,'eivGP')
  if(normalizePath(find.package('eivGP'))!=normalizePath(file.path(p$library,'eivGP')))stop('Wrong eivGP namespace; use a fresh Rscript process')
  manifest<-read.csv(file.path(p$snapshot,'inst/module-manifest.csv'))
  if(!identical(unname(F('eivgp_module_identity')()),manifest$md5))stop('Installed source identity mismatch')
  priors<-p$normalized_priors;priors$dictionary_type<-NULL
  if(!isTRUE(all.equal(F('mixedgp_v030_priors')(priors,p$model_data$p,p$model_data$d,p$model_data$m_vec),p$normalized_priors,tolerance=1e-12)))
    stop('Priors differ after normalization')
  frozen<-readRDS(p$data_source);tr<-frozen$data$train;te<-frozen$data$test
  if(a$action=='validate') {
    cat('Pinned package, normalized priors, frozen source and archived input hashes verified. No fitting.\n')
    return(invisible(p))
  }
  panel<-F('mixedgp_study2_diagnostic_rows')(te$X,te$C,p$panel_points)
  Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  work<-function(method) {
    folder<-file.path(root,method);dir.create(folder,showWarnings=FALSE)
    tryCatch({
      bundle<-file.path(folder,'fit.rds')
      if(file.exists(bundle))saved<-F('mixedgp_read_artifact')(bundle)else {
        attempt<-length(list.files(folder,'^attempt-[0-9]+[.]rds$'))+1L
        atomic(list(started=Sys.time()),file.path(folder,sprintf('attempt-%03d.rds',attempt)))
        options(eivGP.checkpoint_root=file.path(folder,'checkpoints'),eivGP.checkpoint_seconds=900,
          eivGP.checkpoint_iterations=2000,eivGP.checkpoint_after_seconds=0,
          eivGP.checkpoint_script_start=proc.time()[['elapsed']])
        start<-proc.time()[['elapsed']]
        fit<-eivGP::fit_eivgp(tr$X,tr$y,tr$C,U_obs=p$model_data$U_obs,
          latent_dim=p$model_data$d,engine='multivariate',m_vec=p$model_data$m_vec,
          ident=p$model_data$ident,standardize_U=FALSE,kernel=p$kernel$name,matern_nu=p$kernel$matern_nu,
          priors=priors,sampler=method,sampler_control=p$controls[[method]],
          n_iter=p$iterations,burn=p$warmup,n_chains=4L,seed=p$seed,parallel=TRUE,n_cores=4L,
          ess_min=NULL,verbose=FALSE,store_scores=p$store_scores)
        saved<-list(fit=fit,fit_seconds=if(attempt==1L)proc.time()[['elapsed']]-start else NA_real_,attempt=attempt)
        atomic(F('mixedgp_pack_artifact')(saved),bundle)
      }
      fit<-saved$fit
      eq<-function(x,y)isTRUE(all.equal(unname(x),unname(y),check.attributes=FALSE,tolerance=1e-12))
      if(!eq(fit$priors,p$normalized_priors)||!eq(fit$data$X,p$model_data$X)||!eq(fit$data$y,p$model_data$y)||
         !eq(fit$data$C_ord,p$model_data$C_ord)||!eq(fit$data$U_obs,p$model_data$U_obs)||
         fit$data$ident!=p$model_data$ident||fit$control$sampler!=method)stop('Fitted scientific model differs from reference')
      if(!eq(fit$mcmc$chain_initial,p$reference_initial))stop('Initial chain states differ from reference; comparison requires review')
      series<-list(scientific=F('mixedgp_continuation_target_series')(fit,te$X[panel,,drop=FALSE],te$C[panel,,drop=FALSE],NULL,p$n_latent,p$diagnostic_seed,TRUE,4L),
        invariant=F('mixedgp_study2_invariant_series')(fit),raw=F('mixedgp_study2_raw_series')(fit))
      tabs<-lapply(series,F('mixedgp_summarize_diagnostic_series'))
      for(k in names(tabs))write.csv(tabs[[k]],file.path(folder,paste0(k,'-diagnostics.csv')),row.names=FALSE)
      required<-rbind(tabs$scientific,tabs$invariant)
      passes<-vapply(c('scientific','invariant'),function(k)F('mixedgp_diagnostic_table_pass')(tabs[[k]],1.01,400,400,expected_parameters=names(series[[k]])),logical(1))
      precision<-all(is.finite(required$mcse_sd_ratio)&required$mcse_sd_ratio<=.05)
      health<-if(method=='ESS')TRUE else nrow(fit$diagnostics$hmc)==4L&&all(fit$diagnostics$hmc$divergences==0L & is.finite(fit$diagnostics$hmc$step_size) & fit$diagnostics$hmc$step_size>0 & is.finite(fit$diagnostics$hmc$mean_acceptance) & fit$diagnostics$hmc$mean_acceptance>0)
      row<-data.frame(method=method,status=if(all(passes)&&precision&&health)'passed'else'needs_review',
        iterations=p$iterations,warmup=p$warmup,scientific_pass=passes['scientific'],invariant_pass=passes['invariant'],precision_pass=precision,
        worst_rhat=max(required$rhat),min_bulk=min(required$ess_bulk),min_tail=min(required$ess_tail),
        min_scientific_bulk=min(tabs$scientific$ess_bulk),min_invariant_bulk=min(tabs$invariant$ess_bulk),
        fit_seconds=saved$fit_seconds,min_bulk_per_second=min(required$ess_bulk)/saved$fit_seconds,
        hmc_divergences=if(method=='ESS')NA_integer_ else sum(fit$diagnostics$hmc$divergences),
        hmc_mean_acceptance=if(method=='ESS')NA_real_ else mean(fit$diagnostics$hmc$mean_acceptance),attempt=saved$attempt)
      fake<-fit;fake$scientific_ess<-list(target_series=series$scientific,diagnostics=tabs$scientific,scope='full_history')
      F('mixedgp_save_scientific_trace')(fake,file.path(folder,'scientific-trace-sampler'),method)
      if(method!='ESS')write.csv(fit$diagnostics$hmc,file.path(folder,'hmc-diagnostics.csv'),row.names=FALSE)
      loadings<-lapply(fit$mcmc$samples_by_chain$A,function(a)a[,1,1])
      atomic(list(summary=row,panel=panel,initial_states_verified=TRUE),file.path(folder,'assessment.rds'))
      write.csv(row,file.path(folder,'summary.csv'),row.names=FALSE)
      atomic(list(A11=loadings,scientific=series$scientific,invariant=series$invariant),file.path(folder,'trace-series.rds'))
      list(ok=TRUE)
    },error=function(e){writeLines(conditionMessage(e),file.path(folder,'error.txt'));list(ok=FALSE,message=conditionMessage(e))})
  }
  message('Matched pilot: ',workers,' fits concurrently, four chains each')
  results<-parallel::mclapply(p$methods,work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  report()
  failed<-!vapply(results,function(x)is.list(x)&&isTRUE(x$ok),logical(1))
  if(any(failed))stop('Some methods failed; inspect method error.txt files. No missing results count as passes.')
  invisible(results)
}
if(sys.nframe()==0L) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1L])
  study2_sampler_pilot(dirname(dirname(normalizePath(script))))
}
