# Standalone R data preparation. jsonlite is used only for metadata.

# 1. Input-balanced, class-stratified five-fold assignment.
ocean_fold_permutations <- function(x) {
  if(length(x)==1L)return(matrix(x,nrow=1L))
  do.call(rbind,lapply(seq_along(x),function(i)
    cbind(x[i],ocean_fold_permutations(x[-i]))))
}

ocean_fold_standardize <- function(x, floor_sd=0) {
  centered<-sweep(x,2L,colMeans(x),"-")
  scale<-pmax(sqrt(colMeans(centered^2)),floor_sd)
  stopifnot(all(is.finite(scale)),all(scale>0))
  sweep(centered,2L,scale,"/")
}

ocean_assign_folds <- function(cohort,rule,seed) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)
  n<-nrow(cohort);per_class<-n/20L
  stopifnot(n %in% c(500L,1000L),!anyDuplicated(cohort$profile),
    all(table(factor(cohort$C,levels=1:4))==n/4L))
  features<-matrix(0,n,9L);groups<-list()
  quintets<-function(ids)unname(split(ids,rep(seq_len(length(ids)/5L),each=5L)))
  for(class in 1:4) {
    ids<-which(cohort$C==class)
    x<-cbind(cohort$U[ids],log(cohort$G2pressure[ids]),cohort$G2salinity[ids])
    z<-ocean_fold_standardize(x)
    a<-cbind(z,z*z,z[,1]*z[,2],z[,1]*z[,3],z[,2]*z[,3])
    features[ids,]<-ocean_fold_standardize(a,1e-9)
    ordered<-ids[order(cohort$U[ids],cohort$source_row[ids])]
    if(rule=="joint_quintile_pressure") {
      oxygen_groups<-split(ordered,rep(1:5,each=length(ordered)/5L))
      for(group in oxygen_groups) {
        pressure_order<-group[order(cohort$G2pressure[group],cohort$source_row[group])]
        groups<-c(groups,quintets(pressure_order))
      }
    } else if(rule=="oxygen_interleaved")groups<-c(groups,quintets(ordered))
    else stop("Unknown fold rule: ",rule)
  }
  stopifnot(all(lengths(groups)==5L))
  assignment<-integer(n);sums<-array(0,c(5L,9L,4L))
  permutations<-ocean_fold_permutations(1:5)
  for(group in groups) {
    class<-cohort$C[group[1]]
    permutations<-permutations[sample.int(nrow(permutations)),,drop=FALSE]
    scores<-vapply(seq_len(nrow(permutations)),function(i)
      sum((sums[,,class]+features[group[permutations[i,]],,drop=FALSE])^2),numeric(1))
    ordered<-group[permutations[which.min(scores),]]
    assignment[ordered]<-1:5
    sums[,,class]<-sums[,,class]+features[ordered,,drop=FALSE]
  }
  initial<-sum((sums/per_class)^2);accepted<-0L
  for(step in seq_len(20000L)) {
    group<-groups[[sample.int(length(groups),1L)]]
    pair<-group[sample.int(5L,2L)];i<-pair[1];j<-pair[2]
    fi<-assignment[i];fj<-assignment[j];class<-cohort$C[i]
    delta<-features[j,]-features[i,];a<-sums[fi,,class];b<-sums[fj,,class]
    if(sum((a+delta)^2)+sum((b-delta)^2)<sum(a^2)+sum(b^2)-1e-12) {
      sums[fi,,class]<-a+delta;sums[fj,,class]<-b-delta
      assignment[i]<-fj;assignment[j]<-fi;accepted<-accepted+1L
    }
  }
  cohort$fold<-assignment
  stopifnot(all(table(cohort$fold,cohort$C)==per_class))
  list(cohort=cohort,design=list(fold_rule=rule,fold_seed=seed,
    initial_design_score=initial,final_design_score=sum((sums/per_class)^2),
    accepted_input_only_swaps=accepted))
}

# 2. Nested training calibration with class/cell coverage and cutpoint anchors.
ocean_fold_calibration <- function(train,seed) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)
  anchors20<-anchors50<-integer()
  for(cut in c(50,100,150))for(side in list(train$U<cut,train$U>=cut)) {
    ids<-which(side)
    near<-head(ids[order(abs(train$U[ids]-cut),train$source_row[ids])],2L)
    stopifnot(length(near)==2L)
    anchors20<-c(anchors20,near[1]);anchors50<-c(anchors50,near)
  }
  fill<-function(initial,target,min_cell) {
    picked<-unique(initial)
    for(class in 1:4) {
      for(cell in sort(unique(train$cell[train$C==class]))) {
        while(sum(train$cell[picked]==cell)<min_cell) {
          ids<-setdiff(which(train$cell==cell),picked)
          stopifnot(length(ids)>0L)
          picked<-c(picked,ids[sample.int(length(ids),1L)])
        }
      }
      while(sum(train$C[picked]==class)<target) {
        ids<-setdiff(which(train$C==class),picked)
        ratio<-vapply(ids,function(i)
          sum(train$cell[picked]==train$cell[i])/sum(train$cell==train$cell[i]),numeric(1))
        chosen<-order(ratio,runif(length(ids)))[1]
        picked<-c(picked,ids[chosen])
      }
      stopifnot(sum(train$C[picked]==class)==target)
    }
    picked
  }
  scale<-nrow(train)/400L
  o20<-fill(anchors20,20L*scale,2L*scale)
  o50<-fill(c(o20,anchors50),50L*scale,6L*scale)
  train$calib_o20<-as.integer(seq_len(nrow(train)) %in% o20)
  train$calib_o50<-as.integer(seq_len(nrow(train)) %in% o50)
  stopifnot(all(train$calib_o20<=train$calib_o50),
    all(tapply(train$calib_o20,train$C,sum)==20L*scale),
    all(tapply(train$calib_o50,train$C,sum)==50L*scale))
  train
}

# 3. Nested class-balanced test-input masks, using the existing R seed rule.
ocean_fold_test_masks <- function(test,fold) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(202609260L+fold)
  o20<-o50<-rep(FALSE,nrow(test))
  for(class in 1:4) {
    ids<-which(test$C==class);ids<-ids[sample.int(length(ids))];n<-length(ids)
    n20<-round(.2*n)
    n50<-floor(.5*n)+as.integer(n%%2L==1L && (class+fold)%%2L==0L)
    o20[head(ids,n20)]<-TRUE;o50[head(ids,n50)]<-TRUE
  }
  data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=o20,
    observed_U_O50=o50,mask_seed=202609260L+fold)
}

prepare_ocean_folds <- function(dataset="ocean1000",output_dir=NULL,repository_dir=getwd()) {
  dataset<-match.arg(dataset,c("ocean500","ocean1000"))
  if(!requireNamespace("jsonlite",quietly=TRUE))stop("Install the R package jsonlite for metadata.")
  root<-normalizePath(repository_dir,mustWork=TRUE)
  input<-file.path(root,"data/ocean_silicate",sub("ocean","n",dataset))
  if(is.null(output_dir))output_dir<-file.path(root,"reproduction/ocean/prepared",dataset)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir<-normalizePath(output_dir,mustWork=TRUE)
  if(identical(normalizePath(input),output_dir))stop("Choose an output directory separate from the source cohort.")
  meta_file<-file.path(input,if(dataset=="ocean1000")"cohort_metadata.json" else "freeze_manifest.json")
  metadata<-jsonlite::read_json(meta_file,simplifyVector=FALSE)
  cohort<-read.csv(file.path(input,"cohort_manifest.csv"),stringsAsFactors=FALSE)
  cohort$fold<-NULL
  cohort$C<-as.integer(findInterval(cohort$U,c(50,100,150))+1L)
  n<-nrow(cohort);seed<-as.integer(metadata$sample_seed)
  stopifnot(n==as.integer(sub("ocean","",dataset)))
  rule<-if(n==500L)"joint_quintile_pressure" else "oxygen_interleaved"
  folded<-ocean_assign_folds(cohort,rule,seed+if(n==500L)100L else 120L)
  cohort<-folded$cohort
  write.csv(cohort,file.path(output_dir,"cohort_manifest.csv"),row.names=FALSE)
  for(f in 1:5) {
    train<-ocean_fold_calibration(cohort[cohort$fold!=f,,drop=FALSE],seed+200L+f)
    test<-cohort[cohort$fold==f,,drop=FALSE]
    masks<-ocean_fold_test_masks(test,f)
    stopifnot(!length(intersect(train$profile,test$profile)),
      sum(masks$observed_U_O20)==.2*nrow(test),sum(masks$observed_U_O50)==.5*nrow(test),
      all(masks$observed_U_O20<=masks$observed_U_O50))
    write.csv(train,file.path(output_dir,sprintf("fold%02d_train.csv",f)),row.names=FALSE)
    write.csv(test,file.path(output_dir,sprintf("fold%02d_test.csv",f)),row.names=FALSE)
    write.csv(masks,file.path(output_dir,sprintf("fold%02d_test_U_masks.csv",f)),row.names=FALSE)
  }
  metadata[names(folded$design)]<-folded$design
  metadata$data_stage<-"five-fold inputs"
  metadata$train<-.8*n;metadata$test<-n/5L
  metadata$calibration_seed_rule<-"sample_seed+200+fold"
  metadata$test_mask_seed_rule<-"202609260+fold"
  metadata$preparation_implementation<-"base-R-v1"
  metadata$R_version<-as.character(getRversion())
  metadata$RNG_kind<-RNGkind()
  metadata$preparation_script_md5<-unname(tools::md5sum(file.path(root,"data/ocean_silicate/prepare_folds.R")))
  jsonlite::write_json(metadata,file.path(output_dir,"freeze_manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
  message("Prepared ",n,"-profile data using R: ",output_dir)
  invisible(output_dir)
}

if(sys.nframe()==0L) {
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  root<-normalizePath(file.path(dirname(entry),"../.."),mustWork=TRUE)
  args<-commandArgs(TRUE)
  prepare_ocean_folds(if(length(args))args[1] else "ocean1000",
    if(length(args)>1L)args[2] else NULL,repository_dir=root)
}
