# Author: Xinye Chen. ADNI screening orchestration on the unchanged public package.
adni_output_check <- function(path,repo) {
  path<-normalizePath(path.expand(path),mustWork=FALSE);repo<-normalizePath(repo)
  if(startsWith(paste0(path,'/'),paste0(repo,'/'))&&
     !startsWith(paste0(path,'/'),paste0(repo,'/results/')))
    stop('Choose ignored results/ or an external output directory; never write to versioned source/data.')
  invisible(path)
}
adni_screen_api <- function(repo) {
  repo<-normalizePath(repo);lib<-Sys.getenv('EIVGP_R_LIBRARY',file.path(repo,'.r-library'))
  if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
  if(!requireNamespace('eivGP',quietly=TRUE)||packageVersion('eivGP')<'0.4.0.9013')stop('Install this checkout of eivGP first.')
  ns<-asNamespace('eivGP');manifest<-read.csv(system.file('module-manifest.csv',package='eivGP'))
  native<-grepl('[.]cpp$|^Makevars',manifest$module)
  paths<-file.path(repo,'eivGP',ifelse(native,'src','R'),manifest$module)
  if(!identical(unname(tools::md5sum(paths)),manifest$md5))stop('Installed package differs from pinned checkout.')
  e<-new.env(parent=ns)
  for(n in c('mixedgp_budget_yield','mixedgp_budget_worker','mixedgp_budget_enter','mixedgp_budget_cores','eivgp_module_identity'))e[[n]]<-get(n,ns)
  e$repository_dir<-repo;e$api<-e;e$queue<-eivGP::ocean_workflow_api()
  e$atomic_json<-function(value,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);p<-tempfile(tmpdir=dirname(path));on.exit(unlink(p));jsonlite::write_json(value,p,auto_unbox=TRUE,pretty=TRUE,null='null',na='null');stopifnot(file.rename(p,path))}
  e$atomic_csv<-get('ocean_atomic_csv',ns);e$atomic_rds<-get('ocean_atomic_rds',ns)
  sys.source(file.path(repo,'paper/applications/adni-data.R'),e)
  sys.source(file.path(repo,'paper/applications/adni-evaluation.R'),e)
  e
}
adni_screen_args <- function(args) {
  p<-list(cores=5L,plan=FALSE,reuse_only=FALSE)
  for(a in args){
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%c('manifest','cores','plan','reuse_only'))stop('Unknown screening option: ',k)
    if(k%in%c('plan','reuse_only')){if(!v%in%c('true','false'))stop(k,' must be true or false');v<-v=='true'}
    if(k=='cores'){n<-suppressWarnings(as.numeric(v));if(!is.finite(n)||n<1||n!=floor(n)||n>.Machine$integer.max)stop('cores must be a positive integer');v<-as.integer(n)}
    if(k%in%names(p)&&!k%in%c('cores','plan','reuse_only'))stop('Repeated option: ',k)
    p[[k]]<-v
  }
  if(is.null(p$manifest))stop('Supply --manifest=/external/SCREENING.json')
  p
}
adni_screen_settings <- function(m,iter,cores) {
  smoke<-isTRUE(m$smoke)
  get<-function(n,default)if(is.null(m[[n]]))default else m[[n]]
  list(profile=get('baseline_profile',if(smoke)'pilot'else'paper'),run_profile=get('profile',if(smoke)'quick'else'pilot'),
    action='fit',folds=as.integer(unlist(m$folds)),cells=unlist(m$cells),
    cores=as.integer(cores),chains=as.integer(get('chains',4L)),iter=as.integer(iter),burn=as.integer(m$burn),
    first_checkpoint=as.integer(min(iter,get('first_checkpoint',if(smoke)iter else m$burn+250L))),
    checkpoint_step=as.integer(get('checkpoint_step',1000L)),
    draws_per_chain=as.integer(get('draws_per_chain',if(smoke)5L else 250L)),protocol=m$trial_id,
    samplers=get('samplers',list(O0='HMC',O20='HMC',O50='ESS')))
}
adni_screen_manifest <- function(path,api) {
  m<-jsonlite::fromJSON(normalizePath(path),simplifyVector=FALSE)
  if(!m$schema%in%c('adni-screening-v2','adni-application-v1')||!length(m$cases))stop('Unsupported ADNI manifest')
  if(m$schema=='adni-screening-v2'&&!isTRUE(m$smoke)&&
    (m$burn!=500L||!identical(as.integer(unlist(m$folds)),1:5)))stop('Production screening requires all five folds and burn500')
  stages<-unlist(m$stages)
  if(!length(stages)||any(!is.finite(stages))||any(stages!=as.integer(stages))||any(diff(stages)<=0)||
    length(m$burn)!=1L||!is.finite(m$burn)||m$burn<0||m$burn!=as.integer(m$burn)||min(stages)<=m$burn)stop('Invalid stage budgets')
  if(!length(m$folds)||any(!unlist(m$folds)%in%1:5)||anyDuplicated(unlist(m$folds)))stop('Invalid folds')
  if(!is.null(m$chains)&&(length(m$chains)!=1L||m$chains<1||m$chains!=as.integer(m$chains)))stop('Invalid chain count')
  if(!all(unlist(m$cells)%in%c('O0','O20','O50')))stop('Invalid calibration cells')
  if(anyDuplicated(vapply(m$cases,`[[`,'','id')))stop('Repeated Plan identity')
  for(k in names(m$frozen_inputs))if(!identical(digest::digest(file=k,algo='sha256'),m$frozen_inputs[[k]]))stop('Frozen input changed: ',k)
  inputs<-list()
  for(z in m$cases){
    if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',z$id))stop('Invalid dataset identity')
    adni_output_check(z$output,api$repository_dir)
    inp<-api$adni_plan_inputs(z$input)
    for(cell in inp$cells)if(cell$protocol$plan_id!=z$id||any(cell$test$calibration!=0L))stop('Dataset identity or test mask differs')
    inputs[[z$id]]<-inp
  }
  if(length(inputs)>1L)for(i in 2:length(inputs))for(key in names(inputs[[1]]$cells)){
    a<-inputs[[1]]$cells[[key]];b<-inputs[[i]]$cells[[key]]
    if(!identical(a$train$RID,b$train$RID)||!identical(a$test$RID,b$test$RID)||
       !identical(a$train[,c('Y','X1','X2','U1','U2','calibration')],b$train[,c('Y','X1','X2','U1','U2','calibration')])||
       !identical(a$available,b$available))stop('Paired people, anchors or masks differ')
  }
  list(manifest=m,inputs=inputs)
}
adni_screen_tasks <- function(m,group,stage='') {
  rows<-list();add<-function(plan,fold,cell,method,kind)rows[[length(rows)+1L]]<<-data.frame(plan=plan,fold=as.integer(fold),cell=cell,method=method,kind=kind,stage=stage)
  for(z in m$cases)for(f in unlist(m$folds)){
    if(group=='eivGP')for(cell in unlist(m$cells))add(z$id,f,cell,'EIV-GP','eiv')
    if(group=='competitors'){
      add(z$id,f,'OI','UC-GP+LVGP+EzGP','baseline_group')
      for(cell in intersect(unlist(m$cells),c('O20','O50')))for(method in c('UC-GP','LVGP','EzGP'))add(z$id,f,cell,method,'cc')
    }
  }
  if(!length(rows))return(data.frame())
  x<-do.call(rbind,rows);x$id<-paste(x$plan,x$stage,x$kind,x$fold,x$cell,gsub('[+]','_',x$method),sep='__');x
}
adni_screen_receipt <- function(out,task,phase)file.path(out,'receipts',paste0(phase,'_',task$id,'.rds'))
adni_screen_verified <- function(path,identity) {
  if(!file.exists(path))return(FALSE)
  r<-readRDS(path);if(!identical(r$identity,identity))stop('Receipt identity differs: ',path)
  if(length(r$files)&&(!all(file.exists(r$files))||!identical(unname(tools::md5sum(r$files)),r$md5)))stop('Completed output changed: ',path)
  TRUE
}
adni_screen_phase <- function(tasks,phase,cores,api,FUN,status_dir,chains=4L) {
  dir.create(status_dir,recursive=TRUE,showWarnings=FALSE)
  status<-tasks;status$state<-'PENDING';status$detail<-'';status$worker_pid<-NA_integer_;status$seconds<-NA_real_
  outer<-api$queue$ocean_stage_workers(cores,if(phase=='fit')'fit'else'analyze',nrow(tasks),chains)
  publish<-function()api$atomic_csv(status,file.path(status_dir,paste0(phase,'_status.csv')))
  publish()
  # Match the public Ocean executor: resident jobs share the full active-worker pool.
  api$mixedgp_budget_yield(function()api$queue$ocean_schedule(seq_len(nrow(tasks)),rep(1L,nrow(tasks)),outer,
    function(i,allocated)api$mixedgp_budget_worker(function(){
      logpath<-file.path(status_dir,paste0(phase,'_',make.names(tasks$id[i]),'.log'));con<-file(logpath,'at')
      old_message<-getConnection(sink.number(type='message'));sink(con);sink(con,type='message')
      on.exit({sink(old_message,type='message');sink();close(con)},add=TRUE);t<-proc.time()[3]
      err<-tryCatch({FUN(tasks[i,,drop=FALSE],cores);NULL},error=function(e)conditionMessage(e))
      if(!is.null(err))cat('FAILED:',err,'\n')
      list(state=if(is.null(err))'DONE'else'FAILED',detail=if(is.null(err))''else err,seconds=as.numeric(proc.time()[3]-t))
    }),
    function(i,n,pid,active){status$state[i]<<-'RUNNING';status$worker_pid[i]<<-pid;publish()},
    function(i,r,active){if(!is.list(r)||is.null(r$state))r<-list(state='FAILED',detail='Worker exited without a result',seconds=NA_real_)
      status$state[i]<<-r$state;status$detail[i]<<-r$detail;status$seconds[i]<<-r$seconds;publish()}))
  status
}
adni_screen_collect <- function(m,api) {
  pred<-conv<-u<-pros<-list()
  for(z in m$cases){
    for(mark in list.files(z$output,'^EVALUATION_VALIDATION[.]json$',recursive=TRUE,full.names=TRUE)){
      dir<-dirname(mark)
      for(pair in list(c('prediction_metrics.csv','pred'),c('convergence.csv','conv'),c('U_metrics.csv','u'),c('PROSPECTIVE_U_metrics.csv','pros'))){
        p<-file.path(dir,pair[1]);if(!file.exists(p))next
        val<-cbind(plan=z$id,read.csv(p));dest<-get(pair[2]);dest[[length(dest)+1L]]<-val;assign(pair[2],dest)
      }
    }
    files<-list.files(z$output,'^metrics[.]csv$',recursive=TRUE,full.names=TRUE)
    for(p in files){x<-read.csv(p);if('stage'%in%names(x)&&all(x$stage=='BASELINE'))pred[[length(pred)+1L]]<-cbind(plan=z$id,x)}
  }
  bind<-function(x)if(length(x))api$queue$ocean_bind(x)else data.frame()
  tables<-list(prediction=bind(pred),convergence=bind(conv),training_U=bind(u),prospective_U=bind(pros))
  dest<-file.path(m$output_root,'reports');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(k in names(tables))if(nrow(tables[[k]]))api$atomic_csv(tables[[k]],file.path(dest,paste0(k,'_fold_metrics.csv')))
  p<-tables$prediction;summary<-list()
  if(nrow(p))for(ix in split(seq_len(nrow(p)),interaction(p$plan,p$stage,p$cell,p$method,p$protocol,drop=TRUE))){
    x<-p[ix,,drop=FALSE];v<-x$RMSE_SD
    summary[[length(summary)+1L]]<-data.frame(x[1,c('plan','stage','cell','method','protocol')],completed_folds=length(unique(x$fold)),
       RMSE_SD_mean=mean(v),RMSE_SD_SD=if(length(v)>1)sd(v)else NA_real_)
  }
  five<-bind(summary);if(nrow(five))api$atomic_csv(five,file.path(dest,'fivefold_prediction.csv'))
  stability<-list()
  for(kind in names(tables)[names(tables)%in%c('prediction','training_U','prospective_U')]){
    tab<-tables[[kind]];if(!nrow(tab)||!'stage'%in%names(tab))next
    stage_names<-paste0('S',unlist(m$stages));if(length(stage_names)<2L)next
    x<-tab[tab$stage==stage_names[1],,drop=FALSE];y<-tab[tab$stage==tail(stage_names,1),,drop=FALSE]
    keys<-if(kind=='prediction')c('plan','fold','cell','method','protocol')else c('plan','fold','cell','coordinate')
    if(nrow(x)&&nrow(y)){
      z<-merge(x[,c(keys,'RMSE_SD'),drop=FALSE],y[,c(keys,'RMSE_SD'),drop=FALSE],by=keys,suffixes=c('_first','_last'))
      z$first_stage<-stage_names[1];z$last_stage<-tail(stage_names,1)
      z$target<-kind;z$delta_last_minus_first<-z$RMSE_SD_last-z$RMSE_SD_first;stability[[kind]]<-z
    }
  }
  stable<-bind(stability);if(nrow(stable))api$atomic_csv(stable,file.path(dest,'stage_stability.csv'))
  lines<-c('# ADNI screening: saved results',paste('Author: Xinye Chen. Trial',m$trial_id),
    'Partial fold counts and CORE failures remain explicit. Budget completion does not establish convergence.',
    if(nrow(five))capture.output(print(five,row.names=FALSE))else'No completed scores yet.')
  writeLines(lines,file.path(dest,'REPORT.md'))
  html<-c('<!doctype html><meta charset="utf-8"><title>ADNI screening</title>',
    '<h1>ADNI screening: saved results</h1><p>Author: Xinye Chen. Partial folds and CORE failures remain explicit; completing a budget does not establish convergence.</p>',
    api$queue$ocean_table_html(five,'Response prediction: fold means and sample SD'),
    api$queue$ocean_table_html(tables$convergence,'Full-chain CORE diagnostics'),
    api$queue$ocean_table_html(stable,'Paired stage stability; descriptive screening only'))
  writeLines(html,file.path(dest,'index.html'));invisible(tables)
}
adni_screen_run <- function(repo,p,group) {
  api<-adni_screen_api(repo);bundle<-adni_screen_manifest(p$manifest,api);m<-bundle$manifest;inputs<-bundle$inputs
  cases<-setNames(m$cases,vapply(m$cases,`[[`,'','id'));api$B<-dirname(normalizePath(p$manifest))
  if(p$plan){
    for(g in if(group=='report')c('eivGP','competitors')else group){x<-adni_screen_tasks(m,g);print(x,row.names=FALSE)}
    cat('Public shared CPU pool:',p$cores,'; fit resident limit:',api$queue$ocean_stage_workers(p$cores,'fit',100L,4L),
      '; full-chain reservations/backfill disabled.\n');return(invisible(bundle))
  }
  root<-m$output_root
  adni_output_check(root,repo)
  dir.create(root,recursive=TRUE,showWarnings=FALSE);lock<-file.path(root,'.adni-application.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('ADNI workflow already active: ',lock)
  writeLines(as.character(Sys.getpid()),file.path(lock,'pid'));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  code<-file.path(repo,'paper/applications',c('adni-workflow.R','adni-data.R','adni-evaluation.R','adni-render.R'))
  base_id<-list(schema='adni-upstream-screening-v2',trial=m$trial_id,modules=api$eivgp_module_identity(),
    application=unname(tools::md5sum(code)),burn=m$burn,stages=unlist(m$stages),smoke=isTRUE(m$smoke),
    chains=if(is.null(m$chains))4L else as.integer(m$chains),
    samplers=if(is.null(m$samplers))list(O0='HMC',O20='HMC',O50='ESS')else m$samplers,
    baseline_profile=if(is.null(m$baseline_profile))if(isTRUE(m$smoke))'pilot'else'paper'else m$baseline_profile,
    draws_per_chain=if(is.null(m$draws_per_chain))if(isTRUE(m$smoke))5L else 250L else m$draws_per_chain)
  identities<-list()
  for(z in m$cases){
    id<-c(base_id,list(plan=z$id,inputs=unname(tools::md5sum(inputs[[z$id]]$files))))
    identities[[z$id]]<-id;path<-file.path(z$output,'workflow_identity.rds')
    if(file.exists(path)){if(!identical(readRDS(path),id))stop('Saved workflow identity differs; use a new output root')}else if(group!='report'&&!p$reuse_only)api$atomic_rds(id,path)
  }
  required<-list()
  for(g in if(group=='report')c('eivGP','competitors')else group){
    for(st in if(g=='eivGP')paste0('S',unlist(m$stages))else'BASELINE'){
      tasks<-adni_screen_tasks(m,g,st)
      for(i in seq_len(nrow(tasks))){t<-tasks[i,];out<-cases[[t$plan]]$output
        required[[length(required)+1L]]<-list(path=adni_screen_receipt(out,t,if(g=='eivGP')'evaluate'else'fit'),identity=identities[[t$plan]])}
    }
  }
  if(group=='report'||p$reuse_only){
    if(!all(vapply(required,function(r)adni_screen_verified(r$path,r$identity),logical(1))))stop('Saved-only prerequisites incomplete; no fitting or continuation started')
    return(adni_screen_collect(m,api))
  }
  api$adni_plan_require();job<-file.path(root,'logs',paste0(group,'_',format(Sys.time(),'%Y%m%d_%H%M%S')))
  old<-options(eivGP.core_events=file.path(job,'core-events'));on.exit(options(old),add=TRUE)
  close<-api$mixedgp_budget_enter(p$cores);on.exit(close(),add=TRUE,after=FALSE);cores<-api$mixedgp_budget_cores()
  writeLines(paste('Effective shared worker budget',cores,'; public resident-fit scheduling.'),file.path(root,paste0(group,'_budget.txt')))
  stages<-if(group=='eivGP')unlist(m$stages)else 0L
  for(iter in stages){
    stage<-if(group=='eivGP')paste0('S',iter)else'BASELINE';tasks<-adni_screen_tasks(m,group,stage)
    work<-function(t,n,evaluate=FALSE){
      out<-cases[[t$plan]]$output;inp<-inputs[[t$plan]];id<-identities[[t$plan]];s<-adni_screen_settings(m,if(iter)iter else max(unlist(m$stages)),n)
      phase<-if(evaluate)'evaluate'else'fit';receipt<-adni_screen_receipt(out,t,phase)
      if(adni_screen_verified(receipt,id))return(invisible(NULL))
      if(t$kind=='eiv'){
        path<-api$adni_plan_fit_path(out,t)
        if(evaluate){
          if(!file.exists(path)||api$adni_plan_read_fit(path)$checkpoint$iteration!=iter)stop('Exact stage checkpoint unavailable for evaluation')
          api$analyze_cell(t,inp,out,s,stage);folder<-file.path(out,stage,'analysis',sprintf('fold%02d',t$fold),t$cell)
        }else{
          if(file.exists(path)){
            fit<-api$adni_plan_read_fit(path)
            if(fit$checkpoint$iteration>iter)stop('Cannot refit an earlier stage after continuation; restore its completed evaluation receipt')
            if(fit$checkpoint$burn!=s$burn||fit$control$sampler!=s$samplers[[t$cell]])stop('Checkpoint sampler/burn differs')
          }
          api$adni_plan_fit_task(t,inp,out,s,extend=FALSE);folder<-NULL
        }
      }else{
        methods<-if(t$kind=='baseline_group')c('UC-GP','LVGP','EzGP')else t$method
        failures<-character()
        for(method in methods){one<-t;one$method<-method;tryCatch(api$baseline_cell(one,inp,out,s,cc=t$kind=='cc'),error=function(e)failures<<-c(failures,paste(method,conditionMessage(e))))}
        if(length(failures))stop(paste(failures,collapse='; '))
        folder<-if(t$kind=='cc')file.path(out,'CC',t$cell,'baselines',sprintf('fold%02d',t$fold),t$method)else file.path(out,'baselines',sprintf('fold%02d',t$fold))
      }
      files<-if(is.null(folder))character()else sort(list.files(folder,recursive=TRUE,full.names=TRUE))
      api$atomic_rds(list(identity=id,task=t,files=files,md5=unname(tools::md5sum(files)),finished=as.character(Sys.time())),receipt)
    }
    nchains<-if(is.null(m$chains))4L else as.integer(m$chains)
    status<-adni_screen_phase(tasks,'fit',cores,api,function(t,n)work(t,n),file.path(job,stage),nchains)
    if(group=='eivGP'&&any(status$state=='FAILED'))stop('EIV stage incomplete; checkpoints retained. Inspect ',job)
    if(group=='eivGP'){
      status<-adni_screen_phase(tasks,'evaluate',cores,api,function(t,n)work(t,n,TRUE),file.path(job,stage),nchains)
      if(any(status$state=='FAILED'))stop('EIV evaluation incomplete; continuation not started. Inspect ',job)
    }else if(any(status$state=='FAILED'))warning('Comparator failures retained; saved results remain available.')
    adni_screen_collect(m,api)
  }
  invisible(adni_screen_collect(m,api))
}
adni_screen_main <- function(group,args=commandArgs(TRUE)) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);repo<-normalizePath(file.path(dirname(script),'../..'))
  if(identical(args,'--help')){cat('Rscript paper/applications/adni-',group,'.R --manifest=/external/SCREENING.json --cores=5 [--plan=true] [--reuse-only=true]\n',sep='');return(invisible(NULL))}
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  adni_screen_run(repo,adni_screen_args(args),group)
}

adni_method_main <- function(group,args=commandArgs(TRUE)) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  repo<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
  sys.source(file.path(repo,'paper/applications/adni.R'),envir=environment(adni_method_main))
  if(identical(args,'--help')) {
    cat('Rscript paper/applications/adni-',group,'.R --cores=4 --profile=pilot --data-root=data/adni/inputs --output=results/adni_pilot\n',
      '  --profile=pilot|paper|quick; --config=FILE; --iter=N --burn=N --chains=N --stages=N,N\n',
      '  --folds=1:5 --cells=O0,O20,O50; --plan=true; --reuse-only=true; --report-source=FILE\n',
      '  pilot:4/1500/500; paper:4/50000/5000; overrides are frozen in the run identity.\n',
      '  Run EIV-GP and competitors sequentially; report never fits or extends a chain.\n',sep='')
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  adni_method_run(repo,adni_cli_params(args),group)
}
