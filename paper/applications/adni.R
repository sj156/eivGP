# Author: Xinye Chen. Application settings and routing; scientific kernels are eivGP APIs.
adni_resolve_path <- function(path,repo) {
  path<-path.expand(path)
  if(!grepl('^(/|[A-Za-z]:[/\\\\])',path))path<-file.path(repo,path)
  base<-path;suffix<-character()
  while(!file.exists(base)&&dirname(base)!=base) {
    suffix<-c(basename(base),suffix);base<-dirname(base)
  }
  resolved<-normalizePath(base,mustWork=TRUE)
  if(length(suffix))file.path(resolved,paste(suffix,collapse='/'))else resolved
}
adni_cli_params <- function(args) {
  allowed<-c('cores','profile','dataset','data_root','data_dir','output','output_dir','config','manifest',
    'chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'plan','reuse_only','report_source','report_dir','functional_draws_per_chain','n_latent','run_id')
  p<-list()
  for(a in args) {
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%allowed)stop('Unknown ADNI option: ',k)
    if(k%in%names(p))stop('Repeated option: ',k)
    if(k%in%c('plan','reuse_only','appendix')) {
      if(!v%in%c('true','false'))stop(k,' must be true or false')
      v<-v=='true'
    }
    p[[k]]<-v
  }
  p
}
adni_integer <- function(value,name,minimum=1L) {
  n<-suppressWarnings(as.numeric(value))
  if(length(n)!=1L||!is.finite(n)||n<minimum||n!=floor(n)||n>.Machine$integer.max)
    stop(name,' must be an integer >= ',minimum)
  as.integer(n)
}
adni_integer_list <- function(value,name) {
  if(is.character(value)&&length(value)==1L) {
    if(grepl('^[0-9]+:[0-9]+$',value)) {
      z<-as.integer(strsplit(value,':',fixed=TRUE)[[1]]);value<-seq.int(z[1],z[2])
    }else value<-strsplit(value,',',fixed=TRUE)[[1]]
  }
  unname(vapply(unlist(value),adni_integer,1L,name=name))
}
adni_settings <- function(profile='pilot',cores=4L,...) {
  overrides<-list(...);profile<-match.arg(profile,c('pilot','paper','quick'))
  defaults<-switch(profile,paper=list(chains=4L,iter=50000L,burn=5000L,stages=50000L),
    pilot=list(chains=4L,iter=1500L,burn=500L,stages=c(1000L,1500L)),
    quick=list(chains=1L,iter=80L,burn=30L,stages=80L))
  s<-modifyList(c(defaults,list(profile=profile,cores=cores,folds=if(profile=='quick')1L else 1:5,
    cells=c('O50','O20','O0'),draws_per_chain=if(profile=='quick')10L else 250L,
    checkpoint_step=1000L,baseline_profile=if(profile=='quick')'quick'else'paper',
    sampler_o0='HMC',sampler_o20='HMC',sampler_o50='ESS',representative_fold=1L,
    appendix=TRUE,functional_draws_per_chain=64L,n_latent=256L,run_id='local')),overrides)
  for(n in c('cores','chains','iter','draws_per_chain','checkpoint_step','representative_fold',
             'functional_draws_per_chain','n_latent'))s[[n]]<-adni_integer(s[[n]],n)
  s$burn<-adni_integer(s$burn,'burn',0L);s$folds<-adni_integer_list(s$folds,'folds')
  if(any(!s$folds%in%1:5)||anyDuplicated(s$folds))stop('folds must be distinct members of1:5')
  s$cells<-if(length(s$cells)==1L)strsplit(s$cells,',',fixed=TRUE)[[1]]else unlist(s$cells)
  if(!length(s$cells)||any(!s$cells%in%c('O0','O20','O50'))||anyDuplicated(s$cells))stop('Invalid cells')
  if(!is.null(overrides$stages)) {
    s$stages<-adni_integer_list(overrides$stages,'stages')
    if(!is.null(overrides$iter)&&tail(s$stages,1)!=s$iter)stop('iter must equal the last stage')
    s$iter<-tail(s$stages,1)
  }else if(!is.null(overrides$iter))s$stages<-s$iter
  if(any(diff(s$stages)<=0)||min(s$stages)<=s$burn)stop('Stages must increase and exceed burn')
  s$first_checkpoint<-adni_integer(if(is.null(s$first_checkpoint))min(s$stages[1],s$burn+250L)else s$first_checkpoint,'first_checkpoint')
  if(s$first_checkpoint<=s$burn||s$first_checkpoint>min(s$stages))stop('first-checkpoint must exceed burn and not exceed first stage')
  if(!s$baseline_profile%in%c('paper','quick'))stop('baseline-profile must be paper or quick')
  for(n in c('sampler_o0','sampler_o20','sampler_o50'))if(!s[[n]]%in%c('HMC','ESS'))stop('Unknown sampler')
  attr(s,'explicit_selection')<-intersect(names(overrides),c('folds','cells'))
  s
}
adni_load_api <- function(repo) {
  e<-new.env(parent=globalenv());sys.source(file.path(repo,'paper/applications/adni-workflow.R'),e)
  api<-e$adni_screen_api(repo);api$application<-e;api
}
adni_method_run <- function(repo,params,group) {
  if(!is.null(params$config)) {
    config<-jsonlite::fromJSON(adni_resolve_path(params$config,repo),simplifyVector=TRUE)
    params<-modifyList(config,params)
  }
  get<-function(n,default)if(is.null(params[[n]]))default else params[[n]]
  if(!is.null(params$manifest)) {
    env<-new.env(parent=globalenv());sys.source(file.path(repo,'paper/applications/adni-workflow.R'),env)
    return(env$adni_screen_run(repo,list(manifest=adni_resolve_path(params$manifest,repo),
      cores=adni_integer(get('cores',4L),'cores'),plan=isTRUE(params$plan),reuse_only=isTRUE(params$reuse_only)),group))
  }
  fields<-c('chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'functional_draws_per_chain','n_latent','run_id')
  s<-do.call(adni_settings,c(list(profile=get('profile','pilot'),cores=get('cores',4L)),params[intersect(names(params),fields)]))
  dataset<-get('dataset','adni')
  if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',dataset))stop('Invalid dataset name')
  default_input<-if(dataset=='adni')file.path('data','adni','inputs')else file.path('data','adni',dataset,'inputs')
  input<-adni_resolve_path(get('data_root',get('data_dir',default_input)),repo)
  out<-adni_resolve_path(get('output',get('output_dir',file.path('results',paste(dataset,s$profile,sep='_')))),repo)
  if(isTRUE(params$plan)&&!dir.exists(file.path(input,'MC'))) {
    cat('Dataset:',dataset,'; data:',input,'\n');print(s[c('profile','chains','iter','burn','cores')])
    cat('Supply private inputs before fitting. No computation started.\n');return(invisible(s))
  }
  report_source<-if(!is.null(params$report_source))adni_resolve_path(params$report_source,repo)else NULL
  saved_only<-group=='report'||isTRUE(params$reuse_only)
  if(saved_only&&!isTRUE(params$plan))return(adni_render_report(repo,input,out,s,report_source,
    report_dir=adni_resolve_path(get('report_dir',out),repo),group=group))
  api<-adni_load_api(repo);inp<-api$adni_plan_inputs(input);p<-inp$cells[[1]]$protocol
  api$application$adni_output_check(out,repo)
  files<-setNames(lapply(inp$files,function(x)digest::digest(file=x,algo='sha256')),inp$files)
  m<-list(schema='adni-application-v1',author='Xinye Chen',profile=s$profile,
    package_version=as.character(packageVersion('eivGP')),
    baseline_profile=s$baseline_profile,trial_id=dataset,smoke=s$profile=='quick',
    output_root=file.path(out,'workflow'),cores=s$cores,chains=s$chains,stages=s$stages,burn=s$burn,
    first_checkpoint=s$first_checkpoint,checkpoint_step=s$checkpoint_step,draws_per_chain=s$draws_per_chain,
    samplers=list(O0=s$sampler_o0,O20=s$sampler_o20,O50=s$sampler_o50),folds=s$folds,cells=s$cells,
    cases=list(list(id=p$plan_id,input=input,output=out,C=p$C_order,N=inp$N,q=p$q)),frozen_inputs=files)
  if(isTRUE(params$plan)) {
    e<-new.env(parent=globalenv());sys.source(file.path(repo,'paper/applications/adni-workflow.R'),e)
    print(e$adni_screen_tasks(m,group),row.names=FALSE);print(s[c('profile','chains','iter','burn','cores')]);return(invisible(m))
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE);manifest<-file.path(out,'application_manifest.json')
  if(file.exists(manifest)) {
    prior<-jsonlite::fromJSON(manifest,simplifyVector=FALSE);a<-m;b<-prior;a$cores<-b$cores<-NULL
    canonical<-function(x)jsonlite::fromJSON(jsonlite::toJSON(x,auto_unbox=TRUE,null='null'),simplifyVector=FALSE)
    if(!identical(canonical(a),canonical(b)))
      stop('Saved ADNI settings/inputs differ; select a new output directory. ',
        paste(all.equal(canonical(a),canonical(b)),collapse='; '))
  }else api$atomic_json(m,manifest)
  e<-new.env(parent=globalenv());sys.source(file.path(repo,'paper/applications/adni-workflow.R'),e)
  result<-e$adni_screen_run(repo,list(manifest=manifest,cores=s$cores,plan=FALSE,reuse_only=FALSE),group)
  adni_render_report(repo,input,out,s,report_dir=file.path(out,'reports',group),group=group)
  invisible(result)
}
adni_render_report <- function(repo,data_dir,output_dir,settings,report_source=NULL,
                              report_dir=output_dir,group='report',rmd=NULL,filename='adni.html') {
  if(!requireNamespace('rmarkdown',quietly=TRUE)||!rmarkdown::pandoc_available())
    stop('HTML requires rmarkdown and Pandoc on PATH or RSTUDIO_PANDOC')
  env<-new.env(parent=globalenv());sys.source(file.path(repo,'paper/applications/adni-reporting.R'),env)
  analysis<-env$adni_report_prepare(repo,data_dir,output_dir,settings,report_source,report_dir,group)
  render_env<-new.env(parent=globalenv());render_env$.adni_precomputed_analysis<-analysis
  render_env$.adni_repository_dir<-repo
  if(is.null(rmd))rmd<-file.path(repo,'paper/applications/adni.Rmd')
  rmarkdown::render(rmd,output_file=filename,output_dir=report_dir,
    params=list(profile=settings$profile,action='report',cores=settings$cores,data_dir=data_dir,
      output_dir=output_dir,report_source=if(is.null(report_source))''else report_source,
      representative_fold=settings$representative_fold,appendix=settings$appendix),envir=render_env,quiet=TRUE)
  cat('Saved ADNI HTML: ',file.path(report_dir,filename),'\n',sep='')
  invisible(list(html=file.path(report_dir,filename),analysis=analysis))
}
if(sys.nframe()==0L) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  repo<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
  source(file.path(repo,'paper/applications/adni-workflow.R'));adni_method_main('report')
}
