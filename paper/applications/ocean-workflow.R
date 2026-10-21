# Shared CLI orchestration only. All fitting, scoring and plotting use eivGP.
ocean_method_params <- function(args,group) {
  special<-c('reuse-only','plan');flags<-setNames(as.list(rep(FALSE,2)),special)
  for(key in special) {
    ix<-grep(paste0('^--',key,'='),args)
    if(length(ix)>1L)stop('Repeated option: ',key)
    if(length(ix)) {
      value<-sub('^[^=]+=','',args[ix])
      if(!value%in%c('true','false'))stop('--',key,' must be true or false')
      flags[[key]]<-value=='true';args<-args[-ix]
    }
  }
  if(any(grepl('^--(action|modules|import-dirs)=',args)))
    stop('Choose the method script; use --reuse-only=true for saved results or --plan=true to inspect tasks.')
  p<-ocean_cli_params(args)
  p$action<-if(flags$plan)'plan' else if(group=='report'||flags[['reuse-only']])'report' else 'fit'
  p$modules<-switch(group,eivGP=c('eiv','prediction','diagnostics'),
    competitors=c('original_baselines','complete_case_baselines','prediction'),
    report=c('eiv','original_baselines','complete_case_baselines','prediction','diagnostics','surface_comparison'))
  if(is.null(p$run_id))p$run_id<-paste0('ocean_',group,if(p$action=='report')'_saved' else '_run')
  p
}

# Reporting historical fits is permitted; changing their fitting identity is not.
ocean_saved_method_check <- function(output_dir,data_dir,settings) {
  F<-function(n)getFromNamespace(n,'eivGP')
  selected<-F('ocean_select_tasks')(settings)
  required<-unlist(lapply(seq_len(nrow(selected)),function(i)
    F('ocean_path')(output_dir,F('ocean_task_files')(selected[i,,drop=FALSE],settings$profile,'analyze'))))
  missing<-required[!file.exists(required)]
  if(length(missing))stop('Saved-only workflow requires completed fits and predictions; no fitting started. Missing: ',paste(head(missing,8),collapse=', '))
  inputs<-list.files(data_dir,'[.]csv$',full.names=TRUE)
  saved<-F('ocean_path')(output_dir,'data',basename(inputs))
  if(!all(file.exists(saved))||!identical(unname(tools::md5sum(inputs)),unname(tools::md5sum(saved))))
    stop('Prepared inputs differ from the saved run; no fitting started.')
  invisible(selected)
}

ocean_method_run <- function(repository_dir,params,group) {
  F<-function(n)getFromNamespace(n,'eivGP')
  api<-new.env(parent=environment());ocean_load_api(repository_dir,api)
  getp<-function(n,default)if(is.null(params[[n]]))default else params[[n]]
  profile<-getp('profile','paper');dataset<-getp('dataset','ocean500')
  environment(ocean_dataset_dir)<-api;environment(ocean_output_dir)<-api
  data_dir<-ocean_dataset_dir(repository_dir,dataset,getp('data_dir',''))
  out<-ocean_output_dir(repository_dir,dataset,profile,getp('output_dir',''))
  params$profile<-profile;params<-ocean_restore_fit_settings(params,out)
  settings<-do.call(api$ocean_settings,params[intersect(names(params),names(formals(api$ocean_settings)))])
  if(settings$action=='plan') {
    print(api$ocean_select_tasks(settings),row.names=FALSE)
    cat('Group:',group,'; shared core budget:',settings$cores,'\n');return(invisible(settings))
  }
  if(settings$action=='report')ocean_saved_method_check(out,data_dir,settings)
  if(group=='report')return(ocean_run_cli(repository_dir,params[setdiff(names(params),c('chains','iter','burn'))]))
  if(settings$action=='fit') {
    identity<-api$ocean_identity(repository_dir,data_dir,settings)
    idfile<-F('ocean_path')(out,'workflow_identity.rds')
    # Never silently split the two method groups into different output roots.
    if(file.exists(idfile))api$ocean_assert_identity(readRDS(idfile),identity,out)
    api$ocean_storage_init(out)
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(out,'.application.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('Ocean application already running: ',lock)
  writeLines(as.character(Sys.getpid()),file.path(lock,'pid'));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  job<-api$ocean_job_dir(out,settings$run_id)
  old<-options(eivGP.core_events=file.path(job,'core-events'));on.exit(options(old),add=TRUE)
  release<-F('mixedgp_budget_enter')(settings$cores);on.exit(release(),add=TRUE,after=FALSE)
  result<-api$ocean_logged(file.path(job,'execution.log'),{
    if(settings$action=='fit') {
      api$ocean_execute(repository_dir,data_dir,out,settings)
      settings$action<-if(group=='eivGP')'analyze' else 'report'
      settings$run_id<-paste0(settings$run_id,'_scores')
    }
    a<-api$ocean_execute(repository_dir,data_dir,out,settings)
    dest<-F('ocean_path')(out,'reports',group)
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    keep<-function(d) {
      if(!nrow(d))return(d)
      d[(d$variant=='EIV')==(group=='eivGP'),,drop=FALSE]
    }
    prediction<-a$prediction
    prediction$fold_metrics<-keep(prediction$fold_metrics)
    prediction$fold_metrics_SD_Y<-keep(prediction$fold_metrics_SD_Y)
    prediction$summary<-prediction$summary[(prediction$summary$variant=='EIV')==(group=='eivGP'),,drop=FALSE]
    F('ocean_export_prediction_report')(prediction,dest)
    sections<-F('ocean_table_html')(prediction$summary,paste(group,'predictive performance'))
    if(group=='competitors') {
      cc<-F('ocean_cc_prediction_tables')(out,a$inputs,a$prediction,
        cc_output_dir=if(nzchar(getp('cc_output_dir','')))api$ocean_resolve_path(params$cc_output_dir,repository_dir) else out,
        profile=profile,dest=file.path(dest,'complete_case'),expected_folds=settings$folds)
      F('ocean_atomic_csv')(cc$cc,file.path(dest,'complete_case','competitor_fold_metrics.csv'))
      sections<-c(sections,F('ocean_table_html')(cc$cc,'Complete-case competitors: fold scores'))
    } else if(profile!='quick' && all(a$prediction$availability$available)) {
      diagnostics_dir<-if(nzchar(getp('diagnostics_dir','')))api$ocean_resolve_path(params$diagnostics_dir,repository_dir) else F('ocean_path')(out,'diagnostics')
      d<-F('ocean_run_diagnostics')(out,data_dir,report_dir=diagnostics_dir,
        representative_fold=getp('representative_fold',1L),build_prediction=FALSE,
        build_surfaces=FALSE,prediction_report=prediction)
      targets<-F('ocean_build_target_traces')(out,a$inputs,settings)
      F('mixedgp_save_plot_jobs')(list(fig_ocean500_target_traces=F('ocean_target_trace_plot')(targets$traces)),dest,function(n)c(7.1,6.5))
      F('ocean_atomic_csv')(targets$summary,file.path(dest,'target_diagnostics.csv'))
      sections<-c(sections,F('ocean_table_html')(d$diagnostics$u$summary,'Masked oxygen'),
        '<h2>Scientific quantities: chain traces</h2><img src="fig_ocean500_target_traces.png" style="max-width:100%">',
        F('ocean_table_html')(targets$summary,'Target-specific diagnostics'))
    }
    F('ocean_atomic_lines')(c('<!doctype html><meta charset="utf-8">',paste0('<title>Ocean ',group,'</title>'),
      paste0('<h1>Ocean ',group,' — ',dataset,' / ',profile,'</h1>'),
      '<p>Fold means and SDs describe this fixed partition. Workflow completion does not establish convergence.</p>',sections),file.path(dest,'index.html'))
    F('ocean_atomic_rds')(list(group=group,settings=settings,no_refit=params$action=='report',
      input_hashes=tools::md5sum(list.files(data_dir,'[.]csv$',full.names=TRUE)),
      fit_identity=readRDS(F('ocean_path')(out,'workflow_identity.rds')),
      analysis_code=F('eivgp_module_identity')(),time=Sys.time()),file.path(dest,'provenance.rds'))
    list(path=dest,html=file.path(dest,'index.html'))
  })
  cat('Saved ',group,' analysis: ',result$html,'\n',sep='');invisible(result)
}

ocean_method_main <- function(group,args=commandArgs(TRUE)) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  root<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
  sys.source(file.path(root,'paper/applications/ocean.R'),envir=environment(ocean_method_main))
  if(identical(args,'--help')) {
    cat('Rscript paper/applications/ocean-',group,'.R --cores=16 --profile=paper --dataset=ocean500 --output=/path/to/results\n',
      '  --reuse-only=true: require saved fits and predictions; never fit or continue.\n',
      '  --plan=true: list selected tasks. --folds=1:5 --cells=O0,O20,O50 --diagnostics-dir=...\n',
      '  EIV-GP and competitors each fit and analyze their own group; report only assembles saved results.\n',sep='')
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  ocean_method_run(root,ocean_method_params(args,group),group)
}
