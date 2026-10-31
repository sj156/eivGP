# Complete paper workload per machine. All inference/scoring stays in eivGP.
study2_paper_workflow <- function(lane,args=list(),repo_root=getwd(),
                                 runner=eivGP::run_paper_study) {
  lane<-match.arg(lane,c("eivGP","competitors"))
  if(any(c("phase","component","stage","competitors")%in%names(args)))
    stop("This script runs the complete paper workload; select methods if needed, not phase/component/stage")
  action<-if(is.null(args$action))"plan"else args$action
  if(!action%in%c("plan","data","run","report"))stop("action must be plan, data, run or report")
  if(!is.null(args$profile)&&args$profile!="paper")stop("Use profile=paper")
  available<-if(lane=="eivGP")c("FULL","CUT")else c("UC-GP","LVGP","EzGP","UC-GP-CC","LVGP-CC","EzGP-CC")
  selected<-if(is.null(args$methods))available else strsplit(args$methods,",",fixed=TRUE)[[1]]
  if(!length(selected)||anyDuplicated(selected)||any(!selected%in%available))
    stop("Allowed methods: ",paste(available,collapse=", "))
  # Retention is operational; both lanes use the package's compact CUT summaries.
  retain<-if(is.null(args[["retain-draws"]]))"false"else tolower(as.character(args[["retain-draws"]]))
  if(length(retain)!=1L||is.na(retain)||!retain%in%c("true","false"))stop("retain-draws must be true or false")
  args[["retain-draws"]]<-retain
  args$action<-action;args$profile<-"paper"
  if(is.null(args$output))args$output<-file.path(path.expand("~/eivGP-results"),paste0("study2-",lane))
  original<-if(lane=="eivGP")intersect(selected,"FULL")else intersect(selected,c("UC-GP","LVGP","EzGP"))
  # Supplying archives explicitly requests saved FULL evaluation, not new core fits.
  supplied_full<-any(c("full-fit-map","full-fit-root")%in%names(args))
  run_original<-length(original)>0L && !(lane=="eivGP"&&supplied_full)
  if(run_original&&any(c("reps","calib")%in%names(args)))
    stop("The original paper grid is fixed; reps/calib require selecting only CUT or -CC methods, or supplying saved FULL archives")
  common<-c("action","profile","output","data-root","cores","mixing-plots","checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds")
  result<-list()
  if(run_original) {
    a<-args[intersect(names(args),common)]
    a$component<-if(lane=="eivGP")"core"else"core-competitors"
    a$competitors<-if(lane=="eivGP")"none"else paste(original,collapse=",")
    result$original<-runner("study2",a,repo_root)
  }
  extra<-if(lane=="eivGP")selected else selected[grepl("-CC$",selected)]
  if(length(extra)) {
    a<-args;a$methods<-paste(extra,collapse=",");a$component<-"follow-up"
    if(lane=="eivGP"&&"FULL"%in%extra&&!supplied_full&&action%in%c("run","report")) {
      run_dir<-result$original$core$run_dir
      if(is.null(run_dir)||!dir.exists(run_dir))stop("Original EIV-GP run did not supply its saved-fit directory")
      # Restrict to this run, q=4 and the two required calibration sizes.
      files<-list.files(run_dir,pattern="[.]rds$",recursive=TRUE,full.names=TRUE)
      map<-expand.grid(rep=1:50,calib=c(20L,50L),KEEP.OUT.ATTRS=FALSE)
      map$path<-vapply(seq_len(nrow(map)),function(i) {
        pattern<-sprintf("^recovery_fit_primary_rep%03d_cal%03d_s2v18-q4-[[:xdigit:]]+[.]rds$",map$rep[i],map$calib[i])
        found<-files[grepl(pattern,basename(files))]
        if(length(found)!=1L)stop("Need one matching newly saved FULL fit: ",pattern,"; found ",length(found),". No replacement fit will be launched.")
        normalizePath(found,mustWork=TRUE)
      },character(1))
      mapdir<-file.path(path.expand(args$output),"paper-inputs")
      dir.create(mapdir,recursive=TRUE,showWarnings=FALSE)
      mapfile<-file.path(mapdir,"full-fit-map.csv")
      if(file.exists(mapfile)) {
        if(!identical(read.csv(mapfile,stringsAsFactors=FALSE),map))stop("Saved FULL map changed; use a new output root")
      } else write.csv(map,mapfile,row.names=FALSE)
      a[["full-fit-map"]]<-mapfile
    }
    result$calibrated<-runner("study2",a,repo_root)
  }
  invisible(result)
}
study2_lane_cli <- function(lane,script) {
  env<-new.env(parent=globalenv())
  sys.source(file.path(dirname(normalizePath(script)),"common.R"),envir=env)
  args<-list()
  for(x in commandArgs(TRUE)) {
    if(!grepl("^--[^=]+=.+$",x))stop("Use --name=value")
    k<-sub("^--([^=]+)=.*$","\\1",x)
    if(k%in%names(args))stop("Repeated option: ",k)
    args[[k]]<-sub("^--[^=]+=","",x)
  }
  env$numerical_package_check()
  study2_paper_workflow(lane,args,env$NUMERICAL_REPO_ROOT)
}
