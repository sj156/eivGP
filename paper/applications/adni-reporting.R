# Author: Xinye Chen. Saved-only ADNI report data, metrics and figures.
adni_table_html <- function(d,caption='',note='') {
  if(!is.data.frame(d)||!nrow(d))return('<p>No saved result available for this quantity.</p>')
  for(j in seq_along(d))if(is.numeric(d[[j]]))d[[j]]<-round(d[[j]],4)
  # Reuse Ocean's padding, rules, header shading, escaping and bold-mean markup.
  eivGP::ocean_workflow_api()$ocean_table_html(d,caption,note)
}
adni_bind <- function(parts) {
  parts<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,parts)
  if(!length(parts))return(data.frame())
  cols<-unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-NA;x[,cols,drop=FALSE]}))
}
adni_csv <- function(path,required=TRUE) {
  if(!file.exists(path)){if(required)stop('Missing saved report prerequisite: ',path);return(data.frame())}
  read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
}
adni_summary <- function(d,groups,metrics) {
  if(!nrow(d))return(data.frame())
  split_key<-interaction(d[,groups,drop=FALSE],drop=TRUE,lex.order=TRUE)
  adni_bind(lapply(split(d,split_key),function(x){
    out<-x[1,groups,drop=FALSE];out$folds<-length(unique(x$fold))
    for(k in metrics){v<-as.numeric(x[[k]]);out[[paste0(k,'_mean')]]<-mean(v);out[[paste0(k,'_SD')]]<-if(length(v)>1)sd(v)else NA_real_}
    out
  }))
}
adni_C_dictionary <- function(items,levels,labels=NULL) {
  meanings<-c(MEMORY1='Remembering shopping items without a list',
    MEMORY2='Remembering recent events',MEMORY3='Recalling conversations a few days later',
    MEMORY4='Remembering where objects were placed',MEMORY5='Repeating stories or questions',
    MEMORY6='Remembering the current date or day',MEMORY7='Remembering already having told someone something',
    MEMORY8='Remembering appointments, meetings or engagements',
    DIVATT3='Concentrating without distraction from the environment')
  if(length(labels))for(n in names(labels))meanings[n]<-as.character(labels[[n]])
  text<-unname(meanings[items]);text[is.na(text)]<-'Description not supplied in dataset metadata'
  data.frame(Model_column=paste0('C',seq_along(items)),ECog_field=items,
    Item_meaning=text,Native_grades=vapply(levels,function(n)paste(seq_len(n),collapse='/'),''),check.names=FALSE)
}
adni_saved_source <- function(output_dir,source=NULL) {
  if(!is.null(source)&&nzchar(source)) {
    m<-jsonlite::fromJSON(source,simplifyVector=FALSE)
    if(!identical(m$schema,'adni-saved-report-v1'))stop('Unsupported saved-report source')
    return(m)
  }
  p<-file.path(output_dir,'application_manifest.json')
  if(!file.exists(p))stop('Saved report requires application_manifest.json or --report-source=FILE')
  m<-jsonlite::fromJSON(p,simplifyVector=FALSE)
  list(schema='adni-saved-report-v1',output=output_dir,input=m$cases[[1]]$input,
    stage=paste0('S',tail(unlist(m$stages),1)),fit_package_version=m$package_version,fit_profile=m$profile,
    fitting_manifest=p,input_sha256=m$frozen_inputs)
}
adni_report_prepare <- function(repo,data_dir,output_dir,settings,report_source=NULL,report_dir=output_dir,group='report') {
  for(p in c('jsonlite','digest','ggplot2','patchwork','posterior'))if(!requireNamespace(p,quietly=TRUE))stop('Report dependency: ',p)
  api<-adni_load_api(repo);api$application$adni_output_check(report_dir,repo)
  m<-adni_saved_source(output_dir,report_source)
  fitting<-if(!is.null(m$fitting_manifest))jsonlite::fromJSON(m$fitting_manifest,simplifyVector=FALSE)else list()
  for(n in setdiff(c('folds','cells'),attr(settings,'explicit_selection')))if(length(fitting[[n]]))
    settings[[n]]<-if(n=='folds')as.integer(unlist(fitting[[n]]))else unlist(fitting[[n]],use.names=FALSE)
  for(n in c('draws_per_chain','baseline_profile'))if(!is.null(fitting[[n]]))settings[[n]]<-fitting[[n]]
  if(length(fitting$samplers))for(cell in names(fitting$samplers))
    settings[[paste0('sampler_',tolower(cell))]]<-fitting$samplers[[cell]]
  if(!is.null(m$fit_profile))settings$profile<-m$fit_profile
  if(!is.null(m$input))data_dir<-m$input
  if(!is.null(m$output))output_dir<-m$output
  input<-api$adni_plan_inputs(data_dir)
  if(length(m$input_sha256))for(p in names(m$input_sha256))
    if(digest::digest(file=p,algo='sha256')!=m$input_sha256[[p]])stop('Frozen report input changed: ',p)
  # Historical fitting identities are read, not rewritten to match reporting code.
  identity<-file.path(output_dir,'workflow_identity.rds')
  if(file.exists(identity)) {
    id<-readRDS(identity)
    if(!is.null(id$inputs)&&!identical(unname(tools::md5sum(input$files)),id$inputs))stop('Saved fit input identity differs')
  }
  if(is.null(m$stage)) {
    stages<-list.dirs(output_dir,recursive=FALSE,full.names=FALSE)
    stages<-stages[grepl('^S[0-9]+K?$',stages)]
    amount<-as.numeric(gsub('K','',sub('^S','',stages)))*ifelse(grepl('K$',stages),1000,1)
    if(!length(stages))stop('No saved EIV stages')
    m$stage<-stages[which.max(amount)]
  }
  pred<-conv<-up<-um<-pp<-pm<-parameters<-statuses<-paths<-list()
  final<-m$stage
  selected<-expand.grid(fold=settings$folds,cell=settings$cells,stringsAsFactors=FALSE)
  want_eiv<-group!='competitors';want_baseline<-group!='eivGP'
  for(i in seq_len(nrow(selected))) {
    f<-selected$fold[i];cell<-selected$cell[i];z<-input$cells[[paste(f,cell,sep='/')]]
    dir<-file.path(output_dir,final,'analysis',sprintf('fold%02d',f),cell)
    if(want_eiv) {
      d<-adni_csv(file.path(dir,'prediction_metrics.csv'));d<-d[d$protocol%in%c('OI','CC-known'),,drop=FALSE]
      pred[[length(pred)+1L]]<-d
      cv<-adni_csv(file.path(dir,'convergence.csv'));conv[[length(conv)+1L]]<-cv
      par<-adni_csv(file.path(dir,'diagnostic_parameters.csv'));par$fold<-f;par$cell<-cell
      parameters[[length(parameters)+1L]]<-par
      paths[[paste(f,cell,sep='/')]]<-list(fold=f,cell=cell,dir=dir,
        fit=file.path(output_dir,'fits',sprintf('fold%02d',f),cell,'fit.rds'))
      if(!file.exists(paths[[paste(f,cell,sep='/')]]$fit))stop('Missing saved EIV fit; report never fits')
      if(cell!='O0') {
        u<-adni_csv(file.path(dir,'U_predictions.csv'));u$fold<-f;u$cell<-cell
        up[[length(up)+1L]]<-u;um[[length(um)+1L]]<-adni_csv(file.path(dir,'U_metrics.csv'))
        p<-adni_csv(file.path(dir,'PROSPECTIVE_U_predictions.csv'));p$fold<-f;p$cell<-cell
        pp[[length(pp)+1L]]<-p;pm[[length(pm)+1L]]<-adni_csv(file.path(dir,'PROSPECTIVE_U_metrics.csv'))
      }
    }
  }
  if(want_baseline)for(f in settings$folds)for(cell in c('OI',intersect(settings$cells,c('O20','O50'))))for(method in c('UC-GP','LVGP','EzGP')) {
    dir<-if(cell=='OI')file.path(output_dir,'baselines',sprintf('fold%02d',f),method)else
      file.path(output_dir,'CC',cell,'baselines',sprintf('fold%02d',f),method)
    version<-m$fit_package_version
    if(length(m$baseline_overrides))for(b in m$baseline_overrides)
      if(b$fold==f&&b$cell==cell&&b$method==method){dir<-b$directory;version<-b$package_version}
    valid<-all(file.exists(file.path(dir,c('metrics.csv','fit.rds','prediction_draws.rds'))))
    state<-'MISSING';optimization<-'not_available'
    if(valid) {
      d<-adni_csv(file.path(dir,'metrics.csv'));d$method<-paste0(method,if(cell=='OI')''else'-CC');d$cell<-cell
      d$stage<-'BASELINE';pred[[length(pred)+1L]]<-d
      a<-readRDS(file.path(dir,'fit.rds'));optimization<-a$model$optimization_status
      if(is.null(optimization))optimization<-'NOT_EXPOSED'
      state<-if(grepl('nonconverged',optimization))'OPTIMIZER_FLAGGED'else'AVAILABLE'
      rm(a)
    }
    statuses[[length(statuses)+1L]]<-data.frame(fold=f,cell=cell,method=paste0(method,if(cell=='OI')''else'-CC'),
      state=state,optimization_status=optimization,fit_package_version=if(is.null(version))'recorded in saved provenance'else version,
      source_directory=dir)
  }
  p<-input$cells[[1]]$protocol
  z<-input$cells[['1/O50']];cohort<-rbind(z$train,z$test);cohort<-cohort[order(cohort$RID),]
  metadata_path<-file.path(data_dir,'DATASET.json')
  meta<-if(file.exists(metadata_path))jsonlite::fromJSON(metadata_path,simplifyVector=FALSE)else list()
  meta<-modifyList(meta,if(is.null(m$metadata))list()else m$metadata)
  if(is.null(meta))meta<-list()
  if(is.null(meta$title))meta$title<-'ADNI memory: partially observed brain imaging'
  if(is.null(meta$population))meta$population<-if(is.null(p$population))'Population specified in private preparation records'else p$population
  if(is.null(meta$calibration_design))meta$calibration_design<-if(is.null(p$calibration_design))
    'C1-proportional stratified sampling, fixed seeded within-stratum order; O20 nested in O50'else
    paste(p$calibration_design,'; strata:',paste(unlist(p$stratification_columns),collapse=', '))
  analysis<-list(repo=repo,api=api,source=m,input=input,output_dir=output_dir,report_dir=report_dir,
    settings=settings,group=group,stage=final,metadata=meta,protocol=p,cohort=cohort,
    prediction=adni_bind(pred),convergence=adni_bind(conv),parameters=adni_bind(parameters),
    U_predictions=adni_bind(up),U_metrics=adni_bind(um),prospective_predictions=adni_bind(pp),
    prospective_metrics=adni_bind(pm),method_status=adni_bind(statuses),paths=paths)
  if(nrow(analysis$convergence)) {
    nums<-unique(analysis$convergence[,c('iterations','burn','chains')])
    if(nrow(nums)!=1L)stop('Selected report stage mixes chain budgets')
    analysis$numerics<-nums
    analysis$settings$chains<-as.integer(nums$chains)
    analysis$settings$iter<-as.integer(nums$iterations)
    analysis$settings$burn<-as.integer(nums$burn)
  }else analysis$numerics<-data.frame(iterations=NA,burn=NA,chains=NA)
  if(!is.null(m$fitting_manifest)) {
    stages<-paste0('S',unlist(fitting$stages))
    analysis$settings$stages<-as.integer(unlist(fitting$stages))
  }else {
    stages<-paste0('S',unlist(m$evaluation_iterations))
    if(length(m$evaluation_iterations))analysis$settings$stages<-as.integer(unlist(m$evaluation_iterations))
  }
  # Historical folders may use explicit stage labels rather than S<iteration>.
  if(length(m$stage_directories))stages<-unlist(m$stage_directories)
  stability<-list()
  for(st in stages)for(i in seq_len(nrow(selected))) {
    path<-file.path(output_dir,st,'analysis',sprintf('fold%02d',selected$fold[i]),selected$cell[i],'prediction_metrics.csv')
    d<-adni_csv(path,FALSE)
    if(nrow(d)){d<-d[d$protocol=='OI',];d$stage<-st;stability[[length(stability)+1L]]<-d}
  }
  analysis$stage_scores<-adni_bind(stability)
  adni_report_exports(analysis)
  analysis$plots<-adni_report_plots(analysis)
  if(want_eiv) {
    env<-new.env(parent=environment());sys.source(file.path(repo,'paper/applications/adni-functionals.R'),env)
    analysis$functionals<-env$adni_saved_functionals(analysis)
  }
  analysis$artifact_manifest<-adni_paper_assets(analysis)
  analysis
}
adni_report_exports <- function(a) {
  dest<-file.path(a$report_dir,'report-assets');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(n in c('prediction','convergence','parameters','U_metrics','prospective_metrics','method_status','stage_scores'))
    if(nrow(a[[n]]))write.csv(a[[n]],file.path(dest,paste0(n,'.csv')),row.names=FALSE)
  a$api$atomic_json(list(author='Xinye Chen',saved_only=TRUE,new_MCMC=FALSE,new_optimization=FALSE,
    fitting_source=a$source,stage=a$stage,saved_fitting_numerics=a$numerics,report_settings=a$settings,
    input_md5=as.list(tools::md5sum(a$input$files)),analysis_modules=a$api$eivgp_module_identity(),
    reporting_source=as.list(tools::md5sum(file.path(a$repo,'paper/applications',c('adni.Rmd','adni-reporting.R','adni-functionals.R'))))),
    file.path(dest,'REPORT_PROVENANCE.json'))
}
adni_metric_table <- function(a,protocol='OI',full=FALSE) {
  d<-a$prediction[a$prediction$protocol==protocol,,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tab<-adni_summary(d,c('cell','method'),metrics)
  if(!nrow(tab))return('<p>No available saved scores for this protocol.</p>')
  tab$CORE<-vapply(seq_len(nrow(tab)),function(i){
    if(tab$method[i]!='EIV-GP')return('N/A')
    q<-a$convergence[a$convergence$cell==tab$cell[i],];paste0(sum(q$core_pass),'/',nrow(q))},'')
  method<-tab$method
  for(i in seq_len(nrow(tab))) {
    flagged<-if(method[i]=='EIV-GP')tab$CORE[i]!=paste0(tab$folds[i],'/',tab$folds[i])else
      any(a$method_status$cell==tab$cell[i]&a$method_status$method==method[i]&a$method_status$state=='OPTIMIZER_FLAGGED')
    if(flagged)method[i]<-paste0(method[i],'*')
  }
  result<-data.frame(Method=method,Calibration=tab$cell,
    Test=if(protocol=='OI')'(X,C)'else ifelse(tab$method=='EIV-GP','(X,U)','(X,U,C)'),
    Folds=tab$folds,check.names=FALSE)
  shown<-if(full)metrics else c('RMSE','RMSE_SD','CRPS_SD','Coverage95','Width95_SD')
  labels<-c(RMSE='RMSE',RMSE_SD='RMSE/sY',MAE_SD='MAE/sY',CRPS_SD='CRPS/sY',
    Coverage95='Coverage (%)',Width95_SD='Width/sY',IntervalScore95_SD='Interval score/sY')
  best_mask<-matrix(FALSE,nrow(result),length(shown),dimnames=list(NULL,unname(labels[shown])))
  for(k in shown) {
    scale<-if(k=='Coverage95')100 else 1
    values<-sprintf(if(k=='Coverage95')'%.2f (%.2f)'else'%.4f (%.4f)',
      scale*tab[[paste0(k,'_mean')]],scale*tab[[paste0(k,'_SD')]])
    if(k%in%c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','IntervalScore95_SD')) {
      strata<-if(protocol=='OI')rep('shared',nrow(tab))else tab$cell
      for(g in unique(strata)) {
        ix<-which(strata==g);best<-min(tab[[paste0(k,'_mean')]][ix])
        hit<-ix[tab[[paste0(k,'_mean')]][ix]==best]
        # Unconverged methods remain eligible; the asterisk records their status.
        best_mask[hit,labels[[k]]]<-TRUE
      }
    }
    result[[labels[[k]]]]<-values
  }
  result$CORE<-tab$CORE
  attr(result,'best')<-best_mask
  adni_table_html(result,if(protocol=='OI')'Ordinal-input prediction: fold mean (sample SD)'else
    'Complete-case prediction: fold mean (sample SD)')
}
adni_u_table <- function(d) {
  if(!nrow(d))return(data.frame())
  tab<-adni_summary(d,c('cell','coordinate'),c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
  tab$gain_vs_C1_percent<-vapply(seq_len(nrow(tab)),function(i){
    v<-d[d$cell==tab$cell[i]&d$coordinate==tab$coordinate[i],]
    100*(1-sqrt(mean(v$RMSE^2))/sqrt(mean(v$C1_mean_RMSE^2)))},0)
  out<-tab[,c('cell','coordinate','folds')]
  for(k in c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
    out[[k]]<-sprintf('%.4f (%.4f)',tab[[paste0(k,'_mean')]],tab[[paste0(k,'_SD')]])
  out$gain_vs_C1_percent<-sprintf('%+.2f',tab$gain_vs_C1_percent);out
}
adni_target_table <- function(a,all_folds=FALSE) {
  if(is.null(a$functionals))return(data.frame())
  d<-a$functionals$precision
  if(!all_folds)d<-d[d$fold==a$settings$representative_fold,,drop=FALSE]
  z<-a$input$cells[[1]];target<-rep('m(x,c)',nrow(d))
  for(j in seq_along(z$u)) {
    name<-paste0('Mean hidden ',a$protocol$U_order[j]);at<-d$target==name
    label<-if(grepl('FDG',a$protocol$U_order[j]))'FDG'else if(grepl('HIP',a$protocol$U_order[j]))'HIP1000/ICV'else a$protocol$U_order[j]
    target[at]<-paste0(label,' mean (',ifelse(d$cell[at]=='O0','model','raw'),')')
  }
  d$Target<-target
  d<-d[order(match(d$target,unique(d$target)),match(d$cell,c('O0','O20','O50')),d$fold),,drop=FALSE]
  out<-data.frame(Target=d$Target,Calibration=d$cell,`Mean (SD)`=sprintf('%.4f (%.4f)',d$mean,d$posterior_SD),
    `R-hat`=sprintf('%.3f',d$Rhat),`Bulk ESS`=sprintf('%.1f',d$bulk_ESS),`Tail ESS`=sprintf('%.1f',d$tail_ESS),
    MCSE=ifelse(abs(d$MCSE)<1e-4,formatC(d$MCSE,format='e',digits=2),sprintf('%.4f',d$MCSE)),check.names=FALSE)
  if(all_folds)out<-data.frame(Fold=d$fold,out,check.names=FALSE)
  out
}
adni_integration_table <- function(a,all_folds=FALSE) {
  if(is.null(a$functionals))return(data.frame())
  d<-a$functionals$integration
  if(!all_folds)d<-d[d$fold==a$settings$representative_fold,,drop=FALSE]
  data.frame(Fold=d$fold,Calibration=d$cell,Nodes=d$n_latent,`Check nodes`=d$n_latent_check,
    `Mean at nodes`=sprintf('%.5f',d$mean),`Mean at check nodes`=sprintf('%.5f',d$mean_check),
    Difference=sprintf('%+.5f',d$difference),check.names=FALSE)
}
adni_report_plots <- function(a) {
  library(ggplot2);library(patchwork)
  z<-a$input$cells[['1/O50']];d<-a$cohort;p<-a$protocol;plots<-list()
  long<-adni_bind(lapply(seq_along(z$u),function(j)data.frame(U=d[[z$u[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$U_order[j])))
  plots$cohort<-ggplot(long,aes(U,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed imaging reference',y='ADNI memory score',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  xlong<-adni_bind(lapply(seq_along(z$x),function(j)data.frame(X=d[[z$x[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$X_order[j])))
  plots$covariates<-ggplot(xlong,aes(X,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed covariate',y='ADNI memory score',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  shown<-seq_len(min(2L,length(z$c)))
  cu<-adni_bind(lapply(shown,function(k)adni_bind(lapply(seq_along(z$u),function(j)
    data.frame(grade=factor(d[[z$c[k]]],levels=seq_len(p$m_vec[k])),U=d[[z$u[j]]],
      item=paste0('C',k,': ',p$C_order[k]),coordinate=p$U_order[j])))))
  plots$C_U<-ggplot(cu,aes(grade,U))+geom_boxplot(width=.55,outlier.shape=NA,fill='grey90',colour='grey35')+
    geom_point(aes(colour=grade),position=position_jitter(width=.12,height=0,seed=6203L),alpha=.45,size=.9)+
    facet_wrap(~item+coordinate,ncol=length(z$u),scales='free_y',labeller=label_wrap_gen(36))+
    labs(x='Native ECog grade (higher = greater functional difficulty)',y='Measured imaging reference (raw units)')+
    theme_minimal(base_size=11)+theme(legend.position='none')
  counts<-adni_bind(lapply(seq_along(z$c),function(j){t<-as.data.frame(table(factor(d[[z$c[j]]],levels=seq_len(p$m_vec[j]))));
    data.frame(item=p$C_order[j],grade=t$Var1,N=t$Freq)}))
  plots$C_counts<-ggplot(counts,aes(grade,N,fill=grade))+geom_col()+facet_wrap(~item,ncol=3,scales='free_y')+
    labs(x='Native ECog grade',y='Participants')+theme_bw()+theme(legend.position='none')
  u<-a$U_predictions
  if(nrow(u)) {
    u$role<-ifelse(u$calibration==1,'Observed calibration','Masked training U')
    u$C1<-factor(u$C1)
    draw_u<-function(x,all=FALSE){
      ggplot(x,aes(true_std,mean_std,colour=C1))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
        geom_linerange(data=x[x$calibration==0,],aes(ymin=lo_std,ymax=hi_std),alpha=.15)+
        geom_point(aes(shape=role),size=1,alpha=.6)+scale_shape_manual(values=c('Masked training U'=16,'Observed calibration'=1))+
        facet_grid(if(all)fold~coordinate else cell~coordinate,scales='free')+
        labs(x='Reference U (training SD units; evaluation only)',y='Posterior mean (training SD units)',colour=p$C_order[1],shape=NULL)+theme_bw()
    }
    plots$imputation<-draw_u(u[u$fold==a$settings$representative_fold,])
    for(cell in unique(u$cell))plots[[paste0('imputation_',cell)]]<-draw_u(u[u$cell==cell,],TRUE)
    hidden<-u[u$calibration==0,]
    density<-adni_bind(list(data.frame(hidden[,c('cell','coordinate')],value=hidden$true_std,distribution='Reference'),
      data.frame(hidden[,c('cell','coordinate')],value=hidden$mean_std,distribution='Posterior means')))
    plots$input_distribution<-ggplot(density,aes(value,colour=distribution,fill=distribution))+
      geom_density(alpha=.12)+facet_grid(cell~coordinate,scales='free')+
      labs(x='Training SD units',y='Density',colour=NULL,fill=NULL)+theme_bw()
    class_rows<-adni_bind(lapply(split(hidden,interaction(hidden$fold,hidden$cell,hidden$coordinate,hidden$C1,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],C1=as.character(x$C1[1]),N=nrow(x),
        RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)),Coverage95=mean(x$true_U>=x$q025&x$true_U<=x$q975))))
    plots$class_table<-adni_summary(class_rows,c('cell','coordinate','C1'),c('RMSE_SD','Coverage95'))
    common<-hidden[as.character(hidden$common_hidden)%in%c('TRUE','True','1'),]
    common_rows<-adni_bind(lapply(split(common,interaction(common$fold,common$cell,common$coordinate,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],N=nrow(x),RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)))))
    plots$common_hidden<-adni_summary(common_rows,c('cell','coordinate'),'RMSE_SD')
    plots$support<-adni_bind(lapply(split(hidden,interaction(hidden$cell,hidden$coordinate,drop=TRUE)),function(x)
      data.frame(cell=x$cell[1],coordinate=x$coordinate[1],fold_rows=nrow(x),lower_nonpositive=sum(x$q025<=0),
        mean_nonpositive=sum(x$posterior_mean<=0),interval_outside_C_grade='Not applicable: C is not an imaging bin')))
  }
  if(nrow(u))plots<-modifyList(plots,adni_input_distributions(a))
  plots
}

adni_input_distributions <- function(a) {
  library(ggplot2);parts<-list()
  for(key in names(a$paths)) {
    path<-a$paths[[key]];if(path$cell=='O0'||path$fold!=a$settings$representative_fold)next
    z<-a$input$cells[[key]]
    hidden<-which(a$input$cells[[paste(path$fold,'O50',sep='/')]]$train$calibration==0L)
    saved<-readRDS(file.path(path$dir,'U_draws.rds'))
    if(!identical(saved$scale,'raw'))stop('Raw saved U draws required for physical distribution checks')
    ids<-unique(round(seq(1,dim(saved$draws)[1],length.out=min(128,dim(saved$draws)[1]))))
    for(j in seq_along(z$u)) {
      ref<-z$train[[z$u[j]]];mu<-mean(ref);ss<-sd(ref)
      draw<-matrix(saved$draws[ids,hidden,j],nrow=length(ids))
      tag<-data.frame(cell=path$cell,coordinate=z$protocol$U_order[j])
      parts[[length(parts)+1L]]<-data.frame(tag,C1=rep(as.character(z$train$C1[hidden]),each=length(ids)),
        value=as.vector((draw-mu)/ss),distribution='Posterior mixture')
      parts[[length(parts)+1L]]<-data.frame(tag,C1=as.character(z$train$C1[hidden]),value=(ref[hidden]-mu)/ss,distribution='Reference')
      means<-colMeans(matrix(saved$draws[,hidden,j],ncol=length(hidden)))
      parts[[length(parts)+1L]]<-data.frame(tag,C1=as.character(z$train$C1[hidden]),value=(means-mu)/ss,distribution='Posterior means')
    }
  }
  d<-adni_bind(parts);if(!nrow(d))return(list())
  d$C1<-factor(d$C1);reference<-d[d$distribution=='Reference',];post<-d[d$distribution!='Reference',]
  draw_distribution<-function(class=FALSE) {
    g<-ggplot()+geom_histogram(data=reference,aes(value,after_stat(density)),bins=15,
      fill='grey80',colour='white')+geom_rug(data=reference,aes(value),alpha=.3)+
      geom_density(data=post,aes(value,colour=distribution,linetype=distribution),linewidth=.7)+
      scale_linetype_manual(values=c('Posterior mixture'='solid','Posterior means'='dashed'))+
      theme_minimal(base_size=11)+theme(legend.position='bottom')+
      labs(x='Common-hidden U (training SD units)',y='Density',colour=NULL,linetype=NULL)
    if(class)g+facet_grid(cell+coordinate~C1,scales='free')else g+facet_grid(cell~coordinate,scales='free')
  }
  list(input_distribution=draw_distribution(),class_distribution=draw_distribution(TRUE))
}

adni_paper_assets <- function(a) {
  dest<-file.path(a$report_dir,'report-assets','paper');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  plots<-list(fig_adni_cohort=a$plots$cohort,fig_adni_C_support=a$plots$C_counts,
    fig_adni_Y_X=a$plots$covariates,fig_adni_C_U=a$plots$C_U,
    fig_adni_target_traces=a$functionals$trace_plot,fig_adni_mean_functions=a$functionals$mean_plot,
    fig_adni_response_surfaces=a$functionals$surface_plot,fig_adni_U_imputation=a$plots$imputation,
    fig_adni_U_distributions=a$plots$input_distribution,fig_adni_C1_U_distributions=a$plots$class_distribution)
  manifest<-list()
  for(name in names(plots)) {
    file<-paste0(name,'.pdf');available<-!is.null(plots[[name]])
    if(available)ggplot2::ggsave(file.path(dest,file),plots[[name]],width=12,height=if(grepl('support|traces|C1',name))9 else 6,
      device=grDevices::pdf,limitsize=FALSE)
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'not_applicable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(file.path(dest,file)))else'')
  }
  oi<-a$prediction[a$prediction$protocol=='OI',,drop=FALSE]
  cc<-a$prediction[a$prediction$protocol=='CC-known',,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tables<-list(tab_adni_design=a$input$audit,
    tab_adni_C_definitions=adni_C_dictionary(a$protocol$C_order,a$protocol$m_vec,a$metadata$C_labels),
    tab_adni_prediction_OI=adni_summary(oi,c('cell','method'),metrics),
    tab_adni_prediction_CC=adni_summary(cc,c('cell','method'),metrics),tab_adni_training_U=adni_u_table(a$U_metrics),
    tab_adni_prospective_U=adni_u_table(a$prospective_metrics),tab_adni_convergence=a$convergence)
  for(name in names(tables)) {
    d<-tables[[name]];available<-nrow(d)>0L;file<-paste0(name,'.csv')
    if(available) {
      write.csv(d,file.path(dest,file),row.names=FALSE)
      writeLines(adni_table_html(d,name),file.path(dest,paste0(name,'.html')))
    }
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'unavailable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(file.path(dest,file)))else'')
  }
  out<-adni_bind(manifest);write.csv(out,file.path(dest,'artifact_manifest.csv'),row.names=FALSE);out
}
