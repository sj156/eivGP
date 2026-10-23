# Author: Xinye Chen. Read-only scientific targets and conditional response slices.
adni_trace_diagnostic <- function(values,chain) {
  groups<-split(values,chain);n<-min(lengths(groups));x<-do.call(cbind,lapply(groups,head,n))
  data.frame(states_per_chain=n,chains=ncol(x),mean=mean(values),posterior_SD=sd(values),
    Rhat=posterior::rhat(x),bulk_ESS=posterior::ess_bulk(x),tail_ESS=posterior::ess_tail(x),MCSE=posterior::mcse_mean(x))
}
adni_baseline_mean <- function(api,dir,X,C,mvec,method) {
  a<-readRDS(file.path(dir,'fit.rds'));x<-sweep(sweep(X,2,a$center,'-'),2,a$scale,'/')
  if(method=='UC-GP') {
    nd<-getFromNamespace('make_mixedgp_competitor_data','eivGP')(x,C,mvec)
    mean<-as.numeric(stats::predict(a$model$fit,newdata=nd,type='UK',seCompute=FALSE,lightReturn=TRUE)$mean)
  }else if(method=='LVGP')mean<-as.numeric(LVGP::LVGP_predict(cbind(x,C),a$model$fit,MSE_on=0)$Y_hat)
  else mean<-api$adni_ezgp_predict_covariance_consistent(api$adni_ezgp_mixed_matrix(x,C,mvec),a$model)$Y_hat
  a$y_center+a$y_scale*mean
}
adni_saved_functionals <- function(a) {
  library(ggplot2)
  dest<-file.path(a$report_dir,'report-assets','functionals');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  traces<-raw_traces<-curves<-surfaces<-precision<-integration<-list()
  code<-file.path(a$repo,'paper/applications/adni-functionals.R')
  summarize_draws<-function(x)data.frame(mean=colMeans(x),lo=apply(x,2,quantile,.025),hi=apply(x,2,quantile,.975))
  for(key in names(a$paths)) {
    path<-a$paths[[key]];f<-path$fold;cell<-path$cell;z<-a$input$cells[[key]];p<-z$protocol
    cache<-file.path(dest,paste0('fold',f,'_',cell,'.rds'))
    identity<-list(fit=unname(tools::md5sum(path$fit)),inputs=unname(tools::md5sum(a$input$files)),
      code=unname(tools::md5sum(code)),states=a$settings$functional_draws_per_chain,n_latent=a$settings$n_latent,
      analysis_modules=a$api$eivgp_module_identity())
    if(file.exists(cache)) {
      saved<-readRDS(cache)
      if(!identical(saved$identity,identity))stop('Derived functional cache differs; choose a fresh report directory: ',cache)
      result<-saved$result
    }else {
      cat('Read-only functionals fold',f,cell,'\n');flush.console()
      fit<-a$api$adni_plan_read_fit(path$fit)
      stopifnot(fit$checkpoint$iteration==a$numerics$iterations,fit$checkpoint$burn==a$numerics$burn)
      info<-fit$mcmc$mcmc_draw_info
      ids<-unlist(lapply(split(seq_len(nrow(info)),info$chain),function(v)
        v[unique(round(seq(1,length(v),length.out=min(length(v),a$settings$functional_draws_per_chain))))]),use.names=FALSE)
      chains<-info$chain[ids]
      xref<-vapply(z$x,function(n)median(z$train[[n]]),0)
      patterns<-apply(z$train[,z$c,drop=FALSE],1,paste,collapse='/');mode<-sort(names(table(patterns))[table(patterns)==max(table(patterns))])[1]
      cref<-as.numeric(strsplit(mode,'/',fixed=TRUE)[[1]])
      age<-seq(quantile(z$train[[z$x[1]]],.1),quantile(z$train[[z$x[1]]],.9),length.out=7)
      grid<-expand.grid(x=age,C1=seq_len(p$m_vec[1]))
      xx<-matrix(rep(xref,each=nrow(grid)),ncol=length(xref));xx[,1]<-grid$x
      cc<-matrix(rep(cref,each=nrow(grid)),ncol=length(cref));cc[,1]<-grid$C1;colnames(cc)<-p$C_order
      seed<-as.integer((as.double(p$seeds$posterior)+63001)%%(2^30-1)+1)
      md<-eivGP::predict_eivgp(fit,new_X=xx,new_C=cc,target='mean',draw_ids=ids,n_latent=a$settings$n_latent,
        include_process_uncertainty=TRUE,return_components=TRUE,seed=seed)
      cm<-attr(md,'conditional_means',exact=TRUE);stopifnot(is.matrix(cm))
      cv<-a$convergence[a$convergence$fold==f&a$convergence$cell==cell,]
      method_label<-if(nrow(cv)&&all(cv$core_pass))'EIV-GP'else'EIV-GP*'
      curve<-cbind(data.frame(fold=f,cell=cell,method=method_label,grid),summarize_draws(md))
      target_grade<-min(2L,p$m_vec[1]);target<-which(grid$C1==target_grade)[4L]
      mt<-data.frame(fold=f,cell=cell,target='m(X reference, C profile)',chain=chains,
        iteration=fit$checkpoint$burn+info$draw_within_chain[ids],value=cm[,target])
      pr<-cbind(data.frame(fold=f,cell=cell,target='m(X reference, C profile)',units='ADNI_MEM'),adni_trace_diagnostic(cm[,target],chains))
      # Repeat only the fixed target with twice as many integration nodes.
      md2<-eivGP::predict_eivgp(fit,new_X=xx[target,,drop=FALSE],new_C=cc[target,,drop=FALSE],target='mean',
        draw_ids=ids,n_latent=2L*a$settings$n_latent,include_process_uncertainty=FALSE,seed=seed)
      check<-data.frame(fold=f,cell=cell,n_latent=a$settings$n_latent,n_latent_check=2L*a$settings$n_latent,
        mean=mean(cm[,target]),mean_check=mean(md2),difference=mean(md2)-mean(cm[,target]))
      refcell<-a$input$cells[[paste(f,'O50',sep='/')]]
      hidden<-which(refcell$train$calibration==0L&refcell$train$C1==target_grade)
      if(!length(hidden))hidden<-which(refcell$train$calibration==0L)
      for(j in seq_along(z$u)) {
        v<-rowMeans(matrix(fit$mcmc$samples_U[,hidden,j],nrow=nrow(info)))
        anchored<-isTRUE(fit$data$latent_scale_anchored)
        if(anchored)v<-fit$data$U_center[j]+fit$data$U_scale[j]*v
        name<-paste0('Mean hidden ',p$U_order[j])
        mt<-rbind(mt,data.frame(fold=f,cell=cell,target=name,chain=info$chain,
          iteration=fit$checkpoint$burn+info$draw_within_chain,value=v))
        pr<-rbind(pr,cbind(data.frame(fold=f,cell=cell,target=name,
          units=if(anchored)'physical U'else'unanchored model U'),adni_trace_diagnostic(v,info$chain)))
      }
      surface<-data.frame()
      if(cell!='O0') {
        cal<-z$train[z$train$calibration==1,,drop=FALSE]
        ug<-lapply(z$u,function(n)seq(quantile(cal[[n]],.1),quantile(cal[[n]],.9),length.out=11))
        if(length(z$u)>=2L)uv<-expand.grid(U1=ug[[1]],U2=ug[[2]])else uv<-data.frame(U1=ug[[1]],U2=0)
        uu<-matrix(rep(vapply(z$u,function(n)median(cal[[n]]),0),each=nrow(uv)),ncol=length(z$u))
        uu[,1]<-uv$U1;if(length(z$u)>=2L)uu[,2]<-uv$U2;colnames(uu)<-p$U_order
        x<-matrix(rep(xref,each=nrow(uv)),ncol=length(xref))
        fd<-eivGP::predict_eivgp(fit,new_X=x,new_U=uu,new_U_scale='raw',target='surface',draw_ids=ids,
          include_process_uncertainty=TRUE,seed=seed+3000L)
        surface<-cbind(data.frame(fold=f,cell=cell,uv),summarize_draws(fd))
      }
      raw<-adni_bind(lapply(c('V','sigma2','r'),function(n){
        v<-fit$mcmc[[paste0('samples_',n)]]
        if(is.null(v))return(data.frame())
        data.frame(fold=f,cell=cell,target=n,chain=info$chain,iteration=fit$checkpoint$burn+info$draw_within_chain,value=as.numeric(v))
      }))
      result<-list(curves=curve,surfaces=surface,traces=mt,raw_traces=raw,precision=pr,integration=check,
        reference=list(X=xref,C=cref,target_C1=target_grade,common_hidden_N=length(hidden)))
      saveRDS(list(identity=identity,result=result),cache)
      rm(fit,md,md2);gc(FALSE)
    }
    curves[[key]]<-result$curves;surfaces[[key]]<-result$surfaces;traces[[key]]<-result$traces
    raw_traces[[key]]<-result$raw_traces;precision[[key]]<-result$precision;integration[[key]]<-result$integration
  }
  # Native observable-mean curves on the same fixed training-defined slices.
  native<-list()
  for(f in a$settings$folds) {
    z<-a$input$cells[[paste(f,'O50',sep='/')]];p<-z$protocol
    grid<-curves[[paste(f,a$settings$cells[1],sep='/')]][,c('x','C1')]
    xref<-vapply(z$x,function(n)median(z$train[[n]]),0)
    patterns<-apply(z$train[,z$c,drop=FALSE],1,paste,collapse='/');mode<-sort(names(table(patterns))[table(patterns)==max(table(patterns))])[1]
    cref<-as.numeric(strsplit(mode,'/',fixed=TRUE)[[1]])
    xx<-matrix(rep(xref,each=nrow(grid)),ncol=length(xref));xx[,1]<-grid$x
    cc<-matrix(rep(cref,each=nrow(grid)),ncol=length(cref));cc[,1]<-grid$C1
    for(method in c('UC-GP','LVGP','EzGP')) {
      status<-a$method_status[a$method_status$fold==f&a$method_status$cell=='OI'&a$method_status$method==method,]
      if(!nrow(status)||status$state=='MISSING')next
      value<-adni_baseline_mean(a$api,status$source_directory,xx,cc,p$m_vec,method)
      native[[length(native)+1L]]<-data.frame(fold=f,cell='native OI',
        method=paste0(method,if(status$state=='OPTIMIZER_FLAGGED')'*'else''),grid,mean=value,lo=NA_real_,hi=NA_real_)
    }
  }
  out<-list(curves=adni_bind(c(curves,native)),surfaces=adni_bind(surfaces),traces=adni_bind(traces),
    raw_traces=adni_bind(raw_traces),precision=adni_bind(precision),integration=adni_bind(integration))
  for(n in names(out))if(nrow(out[[n]]))write.csv(out[[n]],file.path(dest,paste0(n,'.csv')),row.names=FALSE)
  reference_fold<-a$settings$representative_fold
  plot_curve<-function(x)ggplot(x,aes(x,mean,colour=interaction(method,cell),group=interaction(method,cell)))+
    geom_ribbon(aes(ymin=lo,ymax=hi,fill=interaction(method,cell)),alpha=.06,colour=NA,na.rm=TRUE)+
    geom_line(linewidth=.65)+facet_wrap(~C1,ncol=2)+theme_bw()+
    labs(x=a$protocol$X_order[1],y='m(x,c): ADNI_MEM',colour='Method / calibration',fill='Method / calibration')
  out$mean_plot<-plot_curve(out$curves[out$curves$fold==reference_fold,])
  out$mean_all_plot<-plot_curve(out$curves)+facet_grid(fold~C1)
  if(nrow(out$surfaces)) {
    out$surface_plot<-ggplot(out$surfaces[out$surfaces$fold==reference_fold,],aes(U1,U2,fill=mean))+
      geom_raster()+geom_contour(aes(z=mean),colour='white',alpha=.6)+facet_wrap(~cell)+
      scale_fill_viridis_c()+theme_bw()+labs(x=a$protocol$U_order[1],y=a$protocol$U_order[min(2,a$protocol$d)],fill='f(x,u)')
    out$surface_all_plot<-ggplot(out$surfaces,aes(U1,U2,fill=mean))+geom_raster()+
      facet_grid(fold~cell)+scale_fill_viridis_c()+theme_bw()+labs(x=a$protocol$U_order[1],y=a$protocol$U_order[min(2,a$protocol$d)],fill='f(x,u)')
  }
  traceplot<-function(d)ggplot(d[d$fold==reference_fold,],aes(iteration,value,colour=factor(chain)))+
    geom_line(alpha=.6,linewidth=.25)+facet_grid(target~cell,scales='free_y')+theme_bw()+
    labs(x='Saved iteration',y=NULL,colour='Chain')
  out$trace_plot<-traceplot(out$traces);out$raw_trace_plot<-traceplot(out$raw_traces)
  out
}
