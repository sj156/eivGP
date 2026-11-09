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
  if(is.null(args$cores))args$cores<-"12"
  budget<-suppressWarnings(as.integer(args$cores))
  if(length(budget)!=1L||is.na(budget)||budget<4L||as.character(budget)!=as.character(args$cores))
    stop("Paper reproduction requires an integer cores budget >=4, retaining four chains per fit")
  args$action<-action;args$profile<-"paper"
  if(is.null(args$output))args$output<-file.path(path.expand("~/eivGP-results"),paste0("study2-",lane))
  original<-if(lane=="eivGP")intersect(selected,"FULL")else intersect(selected,c("UC-GP","LVGP","EzGP"))
  # Supplying archives explicitly requests saved FULL evaluation, not new core fits.
  supplied_full<-any(c("full-fit-map","full-fit-root")%in%names(args))
  run_original<-length(original)>0L && !(lane=="eivGP"&&supplied_full)
  if(run_original&&any(c("reps","calib")%in%names(args)))
    stop("The original paper grid is fixed; reps/calib require selecting only CUT or -CC methods, or supplying saved FULL archives")
  common<-c("action","profile","output","data-root","cores","mixing-plots","checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds")
  scientific<-list(outer="32",latent="128",`inner-draws`="4",iter="100000",
    warmup="20000",`inner-iter`="4000",`inner-warmup`="1000",
    `gp-response-latent`="16",`cc-response-latent`="64")
  allowed<-c(common,"methods","reps","calib","full-fit-map","full-fit-root","retain-draws",names(scientific))
  if(length(setdiff(names(args),allowed)))stop("Unknown options: ",paste(setdiff(names(args),allowed),collapse=", "))
  for(k in intersect(names(args),names(scientific)))
    if(!identical(as.character(args[[k]]),scientific[[k]]))stop("Paper fixes ",k,"=",scientific[[k]],"; use the package API for sensitivity studies")
  result<-list()
  if(run_original) {
    a<-args[intersect(names(args),common)]
    a$component<-if(lane=="eivGP")"core"else"core-competitors"
    a$competitors<-if(lane=="eivGP")"none"else paste(original,collapse=",")
    result$original<-runner("study2",a,repo_root)
  }
  extra<-if(lane=="eivGP")selected else selected[grepl("-CC$",selected)]
  paper_cut<-lane=="eivGP" && "CUT"%in%extra && !any(c("reps","calib")%in%names(args))
  if(paper_cut)extra<-setdiff(extra,"CUT")
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
  if(paper_cut) {
    pairs<-utils::read.csv(file.path(repo_root,"paper/manifests/study2/feedback-pairs.csv"))
    if(!identical(names(pairs),c("calib","rep")) || anyDuplicated(pairs) ||
       !setequal(unique(pairs$calib),c(20L,50L)) ||
       !identical(sort(pairs$rep[pairs$calib==20L]),1:28) ||
       !identical(sort(pairs$rep[pairs$calib==50L]),1:14))stop("Invalid published feedback-pair manifest")
    groups<-split(pairs$rep,pairs$calib)
    result$feedback<-lapply(names(groups),function(cal) {
      ids<-groups[[cal]]
      a<-args;a$methods<-"CUT";a$component<-"follow-up";a$reps<-paste(ids,collapse=",")
      a$calib<-cal
      # This remains appendix evidence. Completed CC and FULL comparisons use
      # all 50 datasets, and are never truncated to the CUT subset.
      runner("study2",a,repo_root)
    })
  }
  invisible(result)
}
study2_lane_cli <- function(lane,script) {
  if(lane=="competitors")Sys.setenv(VECLIB_MAXIMUM_THREADS="1",OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
  env<-new.env(parent=globalenv())
  sys.source(file.path(dirname(normalizePath(script)),"numerical-common.R"),envir=env)

  env$numerical_package_check()
  args<-getFromNamespace('mixedgp_parse_cli_args','eivGP')(commandArgs(TRUE))
  if(lane=="competitors") {
    if(is.null(args[["checkpoint-after-seconds"]]))args[["checkpoint-after-seconds"]]<-"0"
    if(is.null(args[["checkpoint-seconds"]]))args[["checkpoint-seconds"]]<-"600"
  }
  study2_paper_workflow(lane,args,env$NUMERICAL_REPO_ROOT)
}

# A second invocation should report an existing verified local run, without
# entering the fitting driver or changing its locks/checkpoints.
study2_active_competitor_run <- function(root,alive=function(pid)
  getFromNamespace("mixedgp_followup_process_alive","eivGP")(pid)) {
  locks<-list.dirs(root,recursive=TRUE,full.names=TRUE)
  locks<-locks[basename(locks)=="RUN.lock"]
  for(lock in locks) {
    path<-file.path(lock,"owner.rds")
    if(!file.exists(path))next
    bytes<-tryCatch(readBin(path,"raw",n=file.info(path)$size),error=function(e)NULL)
    owner<-tryCatch(readRDS(path),error=function(e)NULL)
    pid<-owner$pid
    if(length(pid)!=1L||!is.numeric(pid)||is.na(pid)||!is.finite(pid)||
       pid<1||pid!=floor(pid)||pid>.Machine$integer.max)next
    if(!identical(owner$host,unname(Sys.info()[["nodename"]])))next
    if(!isTRUE(alive(as.integer(pid))))next
    after<-tryCatch(readBin(path,"raw",n=file.info(path)$size),error=function(e)NULL)
    if(identical(bytes,after))return(list(pid=as.integer(pid),lock=lock))
  }
  NULL
}

# Preserve the failed experiment as a read-only source. Copy verified artifacts
# once; subsequent invocations operate only on the resumed output tree.
study2_resume_competitors <- function(args,repo_root,runner=study2_paper_workflow,
                                     check_runtime=function()invisible(NULL)) {
  action<-if(is.null(args$action))"plan"else args$action
  if(!action%in%c("plan","run","report"))stop("Resume action must be plan, run or report")
  candidate<-file.path(repo_root,"results","study2-competitors")
  if(is.null(args[["source-root"]]))args[["source-root"]]<-if(dir.exists(candidate))candidate else file.path(path.expand("~/eivGP-results"),"study2-competitors")
  source<-normalizePath(path.expand(args[["source-root"]]),mustWork=TRUE)
  if(is.null(args$output))args$output<-paste0(source,"-resumed")
  target<-path.expand(args$output)
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  target<-if(dir.exists(target))normalizePath(target,mustWork=TRUE)else file.path(normalizePath(dirname(target),mustWork=TRUE),basename(target))
  within<-function(x,y)identical(x,y)||startsWith(x,paste0(y,"/"))
  if(within(target,source)||within(source,target))stop("Resume output must be separate from the original experiment")
  args[["source-root"]]<-NULL;args$output<-target
  if(is.null(args$cores))args$cores<-"12"
  if(is.null(args[["checkpoint-after-seconds"]]))args[["checkpoint-after-seconds"]]<-"0"
  if(is.null(args[["checkpoint-seconds"]]))args[["checkpoint-seconds"]]<-"600"
  if(is.null(args$methods))args$methods<-"UC-GP-CC,LVGP-CC,EzGP-CC"
  # The original grid is complete on this machine; CC defaults continue only
  # its interrupted calibrated comparisons. Explicit methods may include it.
  if(is.null(args[["data-root"]])) {
    local_data<-file.path(repo_root,"results","data","study2")
    args[["data-root"]]<-if(dir.exists(local_data))local_data else file.path(path.expand("~/eivGP-results"),"data","study2")
  }
  marker<-file.path(target,"resume-source.rds")
  running_status<-function() {
    active<-study2_active_competitor_run(target)
    if(is.null(active))return(NULL)
    message("Already running (PID ",active$pid,"). The existing experiment continues in the background.\n",
      "Outputs: ",target,"\nNo second run was started. Rerun after it exits to resume unfinished work.")
    list(status="already_running",pid=active$pid,output=target)
  }
  if(action=="run" && file.exists(marker)) {
    if(!identical(readRDS(marker)$source,source))stop("This resume output belongs to another source experiment")
    status<-running_status()
    if(!is.null(status))return(invisible(status))
  }
  # An existing worker uses its already loaded runtime. Check the current
  # checkout only when this invocation would start work or build a report.
  check_runtime()
  if(action=="plan") {
    message("Original results (unchanged): ",source,"\nResume output: ",target)
    return(runner("competitors",args,repo_root))
  }
  if(action=="report" && !file.exists(marker))stop("Run the resume script before reporting its output")
  if(!file.exists(marker)) {
    dirs<-list.dirs(source,recursive=TRUE,full.names=TRUE)
    locks<-dirs[grepl("([.]lock$|^LOCK$)",basename(dirs))]
    for(lock in locks) {
      rds<-file.path(lock,"owner.rds");csv<-file.path(lock,"owner.csv");pidfile<-file.path(lock,"pid")
      owner<-if(file.exists(rds))readRDS(rds)else if(file.exists(csv))read.csv(csv)else
        if(file.exists(pidfile))list(pid=suppressWarnings(as.integer(readLines(pidfile,warn=FALSE))))else NULL
      if(length(owner$pid)!=1L||is.na(owner$pid)||owner$pid<1L)stop("Cannot verify source lock: ",lock)
      if(!is.null(owner$host)&&!identical(as.character(owner$host),unname(Sys.info()[["nodename"]])))stop("Source lock belongs to another host: ",lock)
      alive<-getFromNamespace("mixedgp_followup_process_alive","eivGP")(owner$pid)
      if(!identical(alive,FALSE))stop("Source worker is active or cannot be verified: ",lock)
    }
    files<-list.files(source,all.files=TRUE,no..=TRUE,recursive=TRUE,full.names=TRUE)
    relative<-substring(files,nchar(source)+2L)
    keep<-!grepl("(^|/)([^/]*[.]lock|LOCK)(/|$)",relative)
    files<-files[keep];relative<-relative[keep]
    if(!length(files)||any(nzchar(Sys.readlink(files))))stop("Source must contain ordinary saved files")
    before<-setNames(unname(tools::md5sum(files)),relative)
    dir.create(target,recursive=TRUE,showWarnings=FALSE)
    release<-getFromNamespace("mixedgp_followup_lock","eivGP")(file.path(target,".snapshot.lock"))
    on.exit(release(),add=TRUE)
    # An unmarked nonempty output is never adopted accidentally. A partial copy
    # made by this script can resume, verifying each already copied file.
    pending<-file.path(target,"resume-copy.rds")
    if(!file.exists(pending)) {
      existing<-setdiff(list.files(target,all.files=TRUE,no..=TRUE),".snapshot.lock")
      if(length(existing))stop("Resume output already contains unrelated files: ",target)
      saveRDS(list(source=source,files=before),pending)
    } else if(!identical(readRDS(pending),list(source=source,files=before)))stop("Source changed since the partial resume copy")
    message("Copying ",length(files)," saved files; original results stay unchanged")
    for(i in seq_along(files)) {
      dest<-file.path(target,relative[i]);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
      if(file.exists(dest)) {
        if(unname(tools::md5sum(dest))!=before[i])stop("Existing resume copy differs: ",dest)
      } else if(!file.copy(files[i],dest,copy.mode=TRUE,copy.date=TRUE))stop("Cannot copy saved artifact: ",files[i])
    }
    after<-setNames(unname(tools::md5sum(files)),relative)
    copied<-setNames(unname(tools::md5sum(file.path(target,relative))),relative)
    if(!identical(before,after)||!identical(before,copied))stop("Source/copy verification failed; no simulation started")
    getFromNamespace("mixedgp_atomic_rds","eivGP")(list(source=source,files=before,created=Sys.time()),marker)
    unlink(pending);release()
  } else if(!identical(readRDS(marker)$source,source))stop("This resume output belongs to another source experiment")
  if(action=="run") {
    status<-running_status()
    if(!is.null(status))return(invisible(status))
  }
  message("Continuing saved experiment at ",target,"; completed jobs are reused")
  runner("competitors",args,repo_root)
}

study2_resume_competitors_cli <- function(script) {
  Sys.setenv(VECLIB_MAXIMUM_THREADS="1",OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
  env<-new.env(parent=globalenv())
  sys.source(file.path(dirname(normalizePath(script)),"numerical-common.R"),envir=env)
  args<-list()
  for(x in commandArgs(TRUE)) {
    if(!grepl("^--[^=]+=.+$",x))stop("Use --name=value")
    k<-sub("^--([^=]+)=.*$","\\1",x)
    if(k%in%names(args))stop("Repeated option: ",k)
    args[[k]]<-sub("^--[^=]+=","",x)
  }
  study2_resume_competitors(args,env$NUMERICAL_REPO_ROOT,check_runtime=env$numerical_package_check)
}
