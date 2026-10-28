# Author: Xinye Chen. OAI application adapter for the shared prepared-fold workflow.
.oai_repository <- local({
  files<-as.character(unlist(lapply(sys.frames(),function(frame) {
    c(if(is.character(frame$ofile))frame$ofile,if(is.character(frame$file))frame$file)
  }),use.names=FALSE))
  files<-files[basename(files)=='oai-tools.R']
  path<-if(length(files))tail(files,1L)else sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  normalizePath(file.path(dirname(normalizePath(path)),'..'))
})
.oai_prepared_engine <- new.env(parent=globalenv())
sys.source(file.path(.oai_repository,'tools/adni-tools.R'),.oai_prepared_engine)

oai_labels <- function(fields,labels=NULL) {
  values<-unlist(labels,use.names=TRUE)
  vapply(seq_along(fields),function(i) {
    if(fields[i]%in%names(values))as.character(values[[fields[i]]])
    else if(length(values)>=i)as.character(values[i])else fields[i]
  },character(1),USE.NAMES=FALSE)
}

oai_load_api <- function(repo=.oai_repository) {
  repo<-normalizePath(repo)
  e<-.oai_prepared_engine$adni_load_api(repo)
  example<-file.path(repo,'data/oai/DATASET.example.json')
  e$oai_metadata<-if(file.exists(example))jsonlite::fromJSON(example,simplifyVector=FALSE)else list()
  e$adni_default_input<-function(dataset,repo)base::file.path(repo,'data',if(dataset=='oai')'oai'else paste0('oai_',dataset))
  source_files<-e$adni_source_files
  e$adni_source_files<-function(repo)unique(c(source_files(repo),base::file.path(repo,'tools',
    c('oai-tools.R','oai.R','oai-eivGP.R','oai-competitors.R','oai-report.R','oai-prepare-folds.R')),
    base::file.path(repo,'paper/applications/oai.Rmd')))
  settings<-e$adni_settings
  e$adni_settings<-function(profile='pilot',cores=4L,...) {
    overrides<-list(...)
    if(profile=='pilot') {
      if(!any(c('iter','stages')%in%names(overrides)))overrides$stages<-1500L
      if(is.null(overrides$checkpoint_step))overrides$checkpoint_step<-500L
    }
    do.call(settings,c(list(profile=profile,cores=cores),overrides))
  }
  read_inputs<-e$adni_plan_inputs
  e$adni_plan_inputs<-function(data_dir) {
    x<-read_inputs(data_dir)
    candidates<-c(base::file.path(data_dir,'DATASET.json'),base::file.path(dirname(data_dir),'DATASET.json'))
    meta<-candidates[file.exists(candidates)]
    if(length(meta))e$oai_metadata<-jsonlite::fromJSON(meta[1],simplifyVector=FALSE)
    x
  }
  run<-e$adni_method_run
  e$adni_method_run<-function(repo,params,group) {
    if(is.null(params$dataset))params$dataset<-'oai'
    run(repo,params,group)
  }
  e$oai_mask_notice<-function(a=NULL) {
    meta<-if(is.null(a))e$oai_metadata else a$metadata
    selected<-isTRUE(meta$test_U_selection$outcome_selected)||
      grepl('post.hoc|outcome.selected',paste(unlist(meta$test_U_design),collapse=' '),ignore.case=TRUE)
    if(!selected)return('')
    paste('The supplied test-U masks were selected after inspecting held-out results.',
      'CC comparisons are exploratory and conditional on that selection; they are not unbiased confirmatory results.',
      'The workflow preserves the supplied masks and performs no performance-based mask search.')
  }
  e$adni_C_dictionary<-function(items,levels,labels=NULL) {
    encoding<-e$oai_metadata$grade_encoding
    native<-vapply(seq_along(items),function(i) {
      map<-encoding[[paste0('C',i)]]
      if(length(map)==levels[i])paste(unlist(map,use.names=FALSE),collapse=' / ')else 'Supplied native-grade mapping'
    },character(1))
    data.frame(Model_column=paste0('C',seq_along(items)),OAI_field=items,
      Item_meaning=oai_labels(items,labels),Clinical_grades=native,
      Model_codes=vapply(levels,function(n)paste(seq_len(n),collapse='/'),character(1)),check.names=FALSE)
  }
  save_plot<-e$adni_save_plot
  e$adni_save_plot<-function(filename,plot,...) {
    if(grepl('training_U',filename,fixed=TRUE))plot<-plot+ggplot2::labs(
      colour=oai_labels(unlist(e$oai_metadata$C_order),e$oai_metadata$C_labels)[1])
    save_plot(filename,plot,...)
  }
  plots<-e$adni_report_plots
  e$adni_report_plots<-function(a) {
    display<-a;fold<-as.integer(unlist(a$metadata$imputation_display_fold))[1]
    if(length(fold)&&!is.na(fold)&&fold%in%a$settings$folds)display$settings$representative_fold<-fold
    p<-plots(display)
    y<-if(length(a$metadata$Y_label))as.character(a$metadata$Y_label)else as.character(a$metadata$Y)
    if(!length(y))y<-'Response Y'
    u<-if(a$protocol$d==1L)oai_labels(a$protocol$U_order,a$metadata$U_labels)else 'Measured U reference'
    c1<-oai_labels(a$protocol$C_order,a$metadata$C_labels)[1]
    p$cohort<-p$cohort+ggplot2::labs(x=u,y=y,colour=c1)
    p$covariates<-p$covariates+ggplot2::labs(y=y,colour=c1)
    p$C_U<-p$C_U+ggplot2::labs(x='Ordered model code (see native-grade dictionary)',y=u)
    p$C_counts<-p$C_counts+ggplot2::labs(x='Ordered model code (see native-grade dictionary)')
    p
  }
  functionals<-e$adni_saved_functionals
  e$adni_saved_functionals<-function(a) {
    f<-functionals(a)
    y<-if(length(a$metadata$Y_label))as.character(a$metadata$Y_label)else 'Response Y'
    f$precision$units[f$precision$units=='ADNI_MEM']<-y
    x<-oai_labels(a$protocol$X_order,a$metadata$X_labels)[1]
    f$mean_plot<-f$mean_plot+ggplot2::labs(x=x,y=paste('m(x,c):',y))
    f$mean_all_plot<-f$mean_all_plot+ggplot2::labs(x=x,y=paste('m(x,c):',y))
    if(a$protocol$d==1L&&nrow(f$surfaces)) {
      surface<-function(d)ggplot2::ggplot(d,ggplot2::aes(U1,mean,colour=cell,fill=cell))+
        ggplot2::geom_ribbon(ggplot2::aes(ymin=lo,ymax=hi),alpha=.15,colour=NA)+
        ggplot2::geom_line(linewidth=.7)+ggplot2::theme_bw()+
        ggplot2::labs(x=oai_labels(a$protocol$U_order,a$metadata$U_labels),y=paste('f(x,u):',y),colour='Calibration',fill='Calibration')
      f$surface_plot<-surface(f$surfaces[f$surfaces$fold==a$settings$representative_fold,])+ggplot2::facet_wrap(~cell)
      f$surface_all_plot<-surface(f$surfaces)+ggplot2::facet_grid(fold~cell)
    }
    utils::write.csv(f$precision,e$adni_path(a$report_dir,'report-assets/functionals/precision.csv'),row.names=FALSE)
    f
  }
  exports<-e$adni_report_exports
  e$adni_report_exports<-function(a) {
    exports(a)
    file<-e$adni_path(a$report_dir,'report-assets/REPORT_PROVENANCE.json')
    record<-jsonlite::fromJSON(file,simplifyVector=FALSE)
    record$application<-'OAI';record$test_U_design<-a$metadata$test_U_design
    record$test_U_selection<-a$metadata$test_U_selection
    record$test_U_selection_notice<-e$oai_mask_notice(a)
    record$imputation_display_fold<-if(length(a$metadata$imputation_display_fold)&&
      as.integer(a$metadata$imputation_display_fold)%in%a$settings$folds)as.integer(a$metadata$imputation_display_fold)else a$settings$representative_fold
    record$reporting_source<-as.list(tools::md5sum(e$adni_source_files(a$repo)))
    e$atomic_json(record,file)
  }
  assets<-e$adni_paper_assets
  e$adni_paper_assets<-function(a) {
    m<-assets(a);folder<-e$adni_path(a$report_dir,'report-assets/paper')
    for(i in seq_len(nrow(m))) {
      old<-m$file[i];new<-gsub('_adni_','_oai_',old,fixed=TRUE)
      if(old!=new&&file.exists(base::file.path(folder,old)))stopifnot(file.rename(base::file.path(folder,old),base::file.path(folder,new)))
      oh<-sub('[.]csv$','.html',old);nh<-sub('[.]csv$','.html',new)
      if(oh!=nh&&file.exists(base::file.path(folder,oh)))stopifnot(file.rename(base::file.path(folder,oh),base::file.path(folder,nh)))
      m$file[i]<-new;m$label[i]<-gsub('_adni_','_oai_',m$label[i],fixed=TRUE)
      if(file.exists(base::file.path(folder,new)))m$MD5[i]<-unname(tools::md5sum(base::file.path(folder,new)))
    }
    utils::write.csv(m,base::file.path(folder,'artifact_manifest.csv'),row.names=FALSE);m
  }
  e$adni_render_report<-function(repo,data_dir,output_dir,settings,report_source=NULL,
      report_dir=output_dir,group='report',rmd=NULL,filename='oai.html') {
    if(!requireNamespace('rmarkdown',quietly=TRUE)||!rmarkdown::pandoc_available())stop('HTML requires rmarkdown and Pandoc')
    a<-e$adni_report_prepare(repo,data_dir,output_dir,settings,report_source,report_dir,group)
    env<-new.env(parent=globalenv());env$.oai_precomputed_analysis<-a;env$.oai_repository_dir<-repo
    if(is.null(rmd))rmd<-base::file.path(repo,'paper/applications/oai.Rmd')
    if(filename=='adni.html')filename<-'oai.html'
    rmarkdown::render(rmd,output_file=filename,output_dir=report_dir,intermediates_dir=report_dir,
      params=list(profile=settings$profile,action='report',cores=settings$cores,data_dir=data_dir,output_dir=output_dir,
        report_source=if(is.null(report_source))''else report_source,representative_fold=settings$representative_fold,appendix=settings$appendix),envir=env,quiet=TRUE)
    cat('Saved OAI HTML: ',base::file.path(report_dir,filename),'\n',sep='')
    invisible(list(html=base::file.path(report_dir,filename),analysis=a))
  }
  e$oai_settings<-e$adni_settings;e$oai_prepare_folds<-e$adni_prepare_folds
  e$oai_method_run<-e$adni_method_run;e$oai_default_output<-e$adni_default_output
  e$oai_load_api<-function(repo=NULL)e
  e
}

oai_settings <- function(profile='pilot',cores=4L,...)oai_load_api()$oai_settings(profile,cores,...)
oai_prepare_folds <- function(cohort_file,metadata_file,output_dir)oai_load_api()$oai_prepare_folds(cohort_file,metadata_file,output_dir)
oai_method_main <- function(group,args=commandArgs(TRUE)) {
  if(identical(args,'--help')) {
    cat('Rscript tools/oai-',group,'.R --profile=paper --cores=4\n',
      '  --profile=pilot|paper|quick; --dataset=NAME; --data-root=PATH; --output=PATH\n',
      '  --cores=N is the sole resource control; new fits use min(N,4) chains.\n',
      '  --plan=true; --reuse-only=true; --iter=N --burn=N --stages=N,N\n',
      '  --folds=1:5 --cells=O0,O20,O50; --ess-min=100; --report-source=FILE\n',
      '  Run EIV and competitors sequentially; reporting reads saved results.\n',sep='')
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  api<-oai_load_api(.oai_repository)
  api$oai_method_run(.oai_repository,api$adni_cli_params(args),group)
}
