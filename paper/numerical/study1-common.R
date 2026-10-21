# Complete paper workload per machine. All inference/scoring stays in eivGP.
study1_paper_workflow <- function(lane,args=list(),repo_root=getwd(),
                                 runner=eivGP::run_paper_study) {
  lane<-match.arg(lane,c("eivGP","competitors"))
  if(any(c("phase","component","stage","competitors")%in%names(args)))
    stop("This script runs the complete paper workload; select methods if needed, not phase/component/stage")
  if(any(c("full-fit-map","full-fit-root","inner-iter","inner-warmup","inner-draws")%in%names(args)))stop("Study I has no CUT or FULL rescoring workflow")
  action<-if(is.null(args$action))"plan"else args$action
  if(!action%in%c("plan","data","run","report"))stop("action must be plan, data, run or report")
  if(!is.null(args$profile)&&args$profile!="paper")stop("Use profile=paper")
  available<-if(lane=="eivGP")"FULL"else c("UC-GP","LVGP","EzGP","UC-GP-CC","LVGP-CC","EzGP-CC")
  selected<-if(is.null(args$methods))available else strsplit(args$methods,",",fixed=TRUE)[[1]]
  if(!length(selected)||anyDuplicated(selected)||any(!selected%in%available))
    stop("Allowed methods: ",paste(available,collapse=", "))
  args$action<-action;args$profile<-"paper"
  if(is.null(args$output))args$output<-file.path(path.expand("~/eivGP-results"),paste0("study1-",lane))
  original<-if(lane=="eivGP")intersect(selected,"FULL")else intersect(selected,c("UC-GP","LVGP","EzGP"))
  run_original<-length(original)>0L
  if(run_original&&any(c("reps","calib","eta")%in%names(args)))
    stop("The original paper grid is fixed; reps/calib/eta require selecting only -CC methods")
  common<-c("action","profile","output","data-root","cores","mixing-plots","checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds")
  targeted<-c("methods","reps","calib","eta","outer","latent","iter","warmup","seed","cc-response-latent")
  if(length(setdiff(names(args),c(common,targeted))))stop("Unknown options: ",paste(setdiff(names(args),c(common,targeted)),collapse=", "))
  if(!any(grepl("-CC$",selected))&&length(setdiff(names(args),c(common,"methods"))))
    stop("Original-grid methods use fixed paper settings; calibrated controls require -CC methods")
  result<-list()
  if(run_original) {
    a<-args[intersect(names(args),common)]
    a$component<-if(lane=="eivGP")"core"else"core-competitors"
    a$competitors<-if(lane=="eivGP")"none"else paste(original,collapse=",")
    result$original<-runner("study1",a,repo_root)
  }
  extra<-if(lane=="eivGP")character()else selected[grepl("-CC$",selected)]
  if(length(extra)) {
    a<-args;a$methods<-paste(extra,collapse=",");a$component<-"complete-case"
    result$calibrated<-runner("study1",a,repo_root)
  }
  invisible(result)
}
study1_lane_cli <- function(lane,script) {
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
  study1_paper_workflow(lane,args,env$NUMERICAL_REPO_ROOT)
}
