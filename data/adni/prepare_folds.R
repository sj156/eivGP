# Author: Xinye Chen. Base-R C/X-only folds and nested stratified observation masks.
adni_seed <- function(seed,label) {
  h<-digest::digest(paste(seed,label,sep='/'),algo='sha256',serialize=FALSE)
  as.integer(strtoi(substr(h,1,7),16L)%%(2^30-1)+1)
}
adni_nested_masks <- function(d,strata='C1',seed,fractions=c(.2,.5)) {
  key<-do.call(paste,c(d[,strata,drop=FALSE],sep='/'));groups<-split(seq_len(nrow(d)),key)
  ranked<-lapply(names(groups),function(g){set.seed(adni_seed(seed,g));v<-groups[[g]];v[sample.int(length(v))]})
  names(ranked)<-names(groups);counts<-lengths(groups);previous<-integer(length(groups));out<-list()
  for(fraction in fractions) {
    n<-as.integer(floor(nrow(d)*fraction+.5));ideal<-n*counts/sum(counts)
    quota<-pmax(previous,pmin(counts,floor(ideal)))
    while(sum(quota)<n){candidates<-which(quota<counts);j<-candidates[which.max(ideal[candidates]-quota[candidates])];quota[j]<-quota[j]+1L}
    while(sum(quota)>n){candidates<-which(quota>previous);j<-candidates[which.max(quota[candidates]-ideal[candidates])];quota[j]<-quota[j]-1L}
    selected<-unlist(Map(function(v,k)head(v,k),ranked,quota),use.names=FALSE)
    out[[paste0('O',round(100*fraction))]]<-sort(d$RID[selected]);previous<-quota
  }
  stopifnot(all(out$O20%in%out$O50));out
}
adni_balanced_folds <- function(d,C,X,seed) {
  n<-nrow(d);features<-list()
  for(c in C)for(level in sort(unique(d[[c]]))) {
    v<-as.numeric(d[[c]]==level);features[[length(features)+1L]]<-v/sqrt(max(1,sum(v)))
  }
  for(x in X) {
    v<-d[[x]];features[[length(features)+1L]]<-.15*(v-mean(v))/sd(v)
    edges<-unique(quantile(v,c(.2,.4,.6,.8)));bins<-findInterval(v,edges)
    for(level in unique(bins)){z<-as.numeric(bins==level);features[[length(features)+1L]]<-.5*z/sqrt(sum(z))}
  }
  features<-do.call(cbind,features);cap<-rep(n%/%5,5);if(n%%5)cap[seq_len(n%%5)]<-cap[seq_len(n%%5)]+1L
  target<-outer(cap/n,colSums(features));totals<-matrix(0,5,ncol(features));sizes<-integer(5);fold<-integer(n)
  rarity<-apply(sapply(C,function(c)as.integer(table(d[[c]])[as.character(d[[c]])])),1,min)
  set.seed(seed);tie<-runif(n);order<-order(rarity,tie)
  for(i in order) {
    available<-which(sizes<cap)
    cost<-vapply(available,function(k)sum(2*(totals[k,]-target[k,])*features[i,]+features[i,]^2),0)
    best<-available[which.min(cost)];fold[i]<-best;sizes[best]<-sizes[best]+1L;totals[best,]<-totals[best,]+features[i,]
  }
  fold
}
adni_prepare_folds <- function(cohort_file,metadata_file,output_dir) {
  if(!requireNamespace('jsonlite',quietly=TRUE)||!requireNamespace('digest',quietly=TRUE))stop('Install jsonlite and digest')
  if(dir.exists(file.path(output_dir,'MC')))stop('Prepared inputs already exist; never overwrite or reroll frozen folds')
  meta<-jsonlite::fromJSON(metadata_file);d<-read.csv(cohort_file,stringsAsFactors=FALSE)
  for(name in c('fold_seed','calibration_seed','test_seed','model_seed')) {
    v<-meta[[name]]
    if(length(v)!=1L||!is.finite(v)||v<0||v!=floor(v)||v>2^30-1)stop('Supply a fixed integer seed: ',name)
  }
  C<-paste0('C',seq_along(meta$C_order));X<-paste0('X',seq_along(meta$X_order));U<-paste0('U',seq_along(meta$U_order))
  fields<-c('RID','Y',X,U,C)
  if(!all(fields%in%names(d))||anyNA(d[,fields])||any(!is.finite(as.matrix(d[,setdiff(fields,'RID')]))))stop('Invalid prepared cohort columns')
  if(anyDuplicated(d$RID)||nrow(d)<10)stop('Expected unique participant IDs and sufficient rows')
  if(length(meta$m_vec)!=length(C)||length(U)<1||length(C)<length(U))stop('Invalid coordinate dimensions')
  for(j in seq_along(C))if(any(!d[[C[j]]]%in%seq_len(meta$m_vec[j])))stop('C must contain native integer grades')
  for(x in X)if(length(unique(d[[x]]))<=2||sd(d[[x]])<=0)stop('X must be continuous; fix categorical population attributes before preparing data')
  strata<-meta$stratification_columns;if(is.null(strata))strata<-'C1'
  if(any(!strata%in%C))stop('Stratification uses native C only')
  d<-d[order(d$RID),fields]
  original<-read.csv(cohort_file)
  if('fold'%in%names(original)) {
    d$fold<-original$fold[match(d$RID,original$RID)]
    if(any(!d$fold%in%1:5)||!all(1:5%in%d$fold))stop('Invalid supplied folds')
  }else d$fold<-adni_balanced_folds(d,C,X,meta$fold_seed)
  for(f in 1:5)for(j in seq_along(C))if(!all(seq_len(meta$m_vec[j])%in%d[d$fold!=f,C[j]]))
    stop('Every declared native grade must be represented in each training fold; do not silently merge levels or reroll')
  for(f in 1:5) {
    tr<-d[d$fold!=f,fields];te<-d[d$fold==f,fields];tr$calibration<-te$calibration<-0L
    cal<-adni_nested_masks(tr,strata,adni_seed(meta$calibration_seed,paste0('fold',f)))
    mask<-adni_nested_masks(te,strata,adni_seed(meta$test_seed,paste0('fold',f)))
    for(cell in c('O20','O50'))if(qr(cbind(1,as.matrix(tr[tr$RID%in%cal[[cell]],U,drop=FALSE])))$rank<length(U)+1L)
      stop('The frozen calibration vectors have insufficient affine rank; no outcome-guided mask reroll is performed')
    for(cell in c('O0','O20','O50')) {
      dir<-file.path(output_dir,'MC',sprintf('fold%02d',f),cell);dir.create(file.path(dir,'private'),recursive=TRUE)
      train<-tr;train$calibration<-as.integer(train$RID%in%cal[[cell]])
      write.csv(train,file.path(dir,'private/train.csv'),row.names=FALSE);write.csv(te,file.path(dir,'private/test.csv'),row.names=FALSE)
      seed_names<-c('production','posterior','prediction','mixed_prediction','prospective_U','UC_GP','LVGP','EzGP','EzGP_draws')
      seeds<-setNames(lapply(seed_names,function(n)adni_seed(meta$model_seed,paste(f,if(n%in%c('UC_GP','LVGP','EzGP','EzGP_draws'))'baseline'else cell,n,sep='/'))),seed_names)
      p<-list(author='Xinye Chen',plan_id=meta$dataset_id,fold=f,O=cell,N=nrow(d),training_N=nrow(tr),test_N=nrow(te),
        observed_training_U=sum(train$calibration),p=length(X),q=length(C),d=length(U),X_order=meta$X_order,U_order=meta$U_order,
        C_order=meta$C_order,m_vec=meta$m_vec,C_levels=lapply(meta$m_vec,seq_len),
        ident=if(cell=='O0')'lower_triangular'else'none',standardize_U=cell!='O0',test_U_available_RID=if(cell=='O0')integer()else mask[[cell]],
        seeds=seeds,stratification_columns=strata,calibration_design='Proportional stratified sampling; fixed seeded order within stratum; nested20/50',
        fold_design=if('fold'%in%names(original))'Supplied frozen folds'else'C/X-only balanced folds; no Y or U used',
        train_md5=unname(tools::md5sum(file.path(dir,'private/train.csv'))),test_md5=unname(tools::md5sum(file.path(dir,'private/test.csv'))))
      jsonlite::write_json(p,file.path(dir,'PROTOCOL.json'),auto_unbox=TRUE,pretty=TRUE)
    }
  }
  jsonlite::write_json(meta,file.path(output_dir,'DATASET.json'),auto_unbox=TRUE,pretty=TRUE)
  invisible(d[,c('RID','fold')])
}
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(length(args)!=3L)stop('Usage: Rscript data/adni/prepare_folds.R PRIVATE_COHORT.csv PRIVATE_METADATA.json OUTPUT_INPUT_DIR')
  adni_prepare_folds(args[1],args[2],args[3])
}
