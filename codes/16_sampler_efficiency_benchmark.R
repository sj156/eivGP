############################################################
## Bounded computational benchmark; never publication evidence.
## Compares accepted-cache+dense and accepted-cache+Schur for both studies.
##
## Example (run from any directory):
## EIVGP_BENCHMARK_N=48 EIVGP_BENCHMARK_ITER=300 \
## EIVGP_BENCHMARK_BURN=100 EIVGP_BENCHMARK_CHAINS=4 \
## EIVGP_BENCHMARK_CALIB_COUNTS=0,12 EIVGP_BENCHMARK_KERNELS=se \
## OUT_DIR=/absolute/path/to/benchmark-output \
## Rscript revision/codes/16_sampler_efficiency_benchmark.R
##
## Optional: EIVGP_BENCHMARK_STUDIES=1,2; _PRESET=balanced;
## _PARALLEL=true; _CORES=4; _SEED=20260905; _PATH_TOLERANCE=1e-8.
## Default n=60, iter=600, burn=200, chains=4, calib=0,15, kernel=se.
############################################################

arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
codes_dir <- if(length(arg)) dirname(normalizePath(sub("^--file=","",arg[1L]))) else getwd()
repo_dir <- normalizePath(file.path(codes_dir,".."))
read_int <- function(name,default) {
  raw <- Sys.getenv(name,unset=as.character(default))
  value <- suppressWarnings(as.numeric(raw))
  if(length(value)!=1L || !is.finite(value) || value!=as.integer(value)) {
    stop(name," must be one integer.")
  }
  as.integer(value)
}
read_bool <- function(name,default) {
  value <- tolower(Sys.getenv(name,unset=as.character(default)))
  if(value %in% c("true","1","yes")) return(TRUE)
  if(value %in% c("false","0","no")) return(FALSE)
  stop(name," must be true or false.")
}
read_values <- function(name,default) {
  unique(trimws(strsplit(Sys.getenv(name,unset=default),",",fixed=TRUE)[[1L]]))
}
cfg <- list(
  n=read_int("EIVGP_BENCHMARK_N",60L),
  n_iter=read_int("EIVGP_BENCHMARK_ITER",600L),
  burn=read_int("EIVGP_BENCHMARK_BURN",200L),
  n_chains=read_int("EIVGP_BENCHMARK_CHAINS",4L),
  n_cores=read_int("EIVGP_BENCHMARK_CORES",4L),
  seed=read_int("EIVGP_BENCHMARK_SEED",20260905L),
  parallel=read_bool("EIVGP_BENCHMARK_PARALLEL",TRUE),
  preset=Sys.getenv("EIVGP_BENCHMARK_PRESET",unset="balanced"),
  kernels=read_values("EIVGP_BENCHMARK_KERNELS","se"),
  studies=read_values("EIVGP_BENCHMARK_STUDIES","1,2"),
  calibration=suppressWarnings(as.integer(read_values("EIVGP_BENCHMARK_CALIB_COUNTS","0,15"))),
  path_tolerance=suppressWarnings(as.numeric(Sys.getenv("EIVGP_BENCHMARK_PATH_TOLERANCE",unset="1e-8")))
)
if(cfg$n<8L || cfg$n_iter<8L || cfg$burn<0L || cfg$n_iter-cfg$burn<4L ||
   cfg$n_chains<2L || cfg$n_cores<1L || anyNA(cfg$calibration) ||
   any(cfg$calibration<0L | cfg$calibration>cfg$n) ||
   any(!cfg$kernels %in% c("se","matern")) ||
   any(!cfg$studies %in% c("1","2")) ||
   !cfg$preset %in% c("fast","balanced","thorough") ||
   length(cfg$path_tolerance)!=1L || !is.finite(cfg$path_tolerance) ||
   cfg$path_tolerance<=0) stop("Invalid benchmark configuration.")
if(!requireNamespace("posterior",quietly=TRUE)) stop("Benchmark requires package posterior.")

out_dir <- Sys.getenv("EIVGP_BENCHMARK_OUT_DIR",unset=Sys.getenv("OUT_DIR",unset=""))
if(!nzchar(out_dir)) out_dir <- file.path(repo_dir,"tmp","sampler-efficiency-benchmark")
if(grepl("(^|/)(paper-ordinal|tables|figures|publication-results)(/|$)",out_dir)) {
  stop("Choose a benchmark output directory outside manuscript/publication artifacts.")
}
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
out_dir <- normalizePath(out_dir,mustWork=TRUE)
if(file.exists(file.path(out_dir,"benchmark-manifest.rds"))) {
  stop("Output directory already contains a benchmark manifest; choose a new OUT_DIR.")
}
manifest <- list(config=cfg,created_at=Sys.time(),purpose="computational_pilot_not_publication_evidence",
  source_md5=tools::md5sum(file.path(codes_dir,c("00_study1_functions.R",
    "00_study2_functions.R","simulation_helpers.R"))),session=utils::sessionInfo())
saveRDS(manifest,file.path(out_dir,"benchmark-manifest.rds"))

engines <- lapply(c("1","2"),function(study) {
  env <- new.env(parent=globalenv())
  source(file.path(codes_dir,"00_parallel_utils.R"),local=env)
  source(file.path(codes_dir,paste0("00_study",study,"_functions.R")),local=env)
  source(file.path(codes_dir,"simulation_helpers.R"),local=env)
  env
})
names(engines) <- c("1","2")
datasets <- list()
for(study in cfg$studies) {
  e <- engines[[study]]
  datasets[[study]] <- if(study=="1") {
    ## This small computational check conditions on observing each of four
    ## balanced levels; it is not the publication simulation design.
    e$simulate_1d_data(n=cfg$n,n_test=5L,m=4L,threshold_design="balanced",
      min_class_count=1L,seed=cfg$seed+1L)
  } else e$simulate_study2_data(n=cfg$n,n_test=5L,seed=cfg$seed+2L)
}
set.seed(cfg$seed+301L)
calibration_order <- sample.int(cfg$n)

strict_summary <- function(table) {
  required <- c("rhat","ess_bulk","ess_tail")
  valid <- nrow(table)>0L && all(is.finite(as.matrix(table[,required,drop=FALSE])))
  list(finite=valid,
    max_rhat=if(valid) max(table$rhat) else NA_real_,
    min_bulk=if(valid) min(table$ess_bulk) else NA_real_,
    min_tail=if(valid) min(table$ess_tail) else NA_real_)
}
diagnose <- function(fit,study,seconds) {
  e <- engines[[study]]
  raw_series <- if(study=="1") e$mixedgp_study1_raw_series(fit) else e$mixedgp_study2_raw_series(fit)
  raw <- e$mixedgp_summarize_diagnostic_series(raw_series)
  raw$scope <- "raw_coordinates"
  invariant <- if(study=="2") {
    value <- e$mixedgp_summarize_diagnostic_series(e$mixedgp_study2_invariant_series(fit))
    value$scope <- "measurement_invariants"
    value
  } else raw[FALSE,,drop=FALSE]
  selected <- if(study=="2" && !length(fit$data$calib_idx)) invariant else rbind(raw,invariant)
  selected_stats <- strict_summary(selected)
  raw_stats <- strict_summary(raw)
  detail <- rbind(raw,invariant)
  detail$bulk_ess_per_second <- detail$ess_bulk/seconds
  detail$tail_ess_per_second <- detail$ess_tail/seconds
  list(detail=detail,selected=selected,raw=raw,
    selected_stats=selected_stats,raw_stats=raw_stats,
    raw_pass=e$mixedgp_diagnostic_table_pass(raw),
    selected_pass=e$mixedgp_diagnostic_table_pass(selected),
    scope=if(study=="2" && !length(fit$data$calib_idx)) "measurement_invariants_only" else "all_free_raw_and_measurement_coordinates")
}
telemetry <- function(fit) {
  stats <- fit$mcmc$chain_stats
  names_gp <- names(stats)[grepl("^gp_",names(stats))]
  data.frame(metric=names_gp,value=vapply(stats[names_gp],sum,numeric(1)),
             stringsAsFactors=FALSE)
}
finite_fit <- function(fit,study) {
  variables <- if(study=="1") c("samples_u","samples_tau","samples_logtheta","samples_sigma2") else {
    c("samples_U","samples_A","samples_tau","samples_logtheta","samples_sigma2")
  }
  all(vapply(fit$mcmc[variables],function(z) length(z)>0L && all(is.finite(z)),logical(1)))
}
path_difference <- function(a,b,study) {
  names <- if(study=="1") c("samples_u","samples_tau","samples_logtheta","samples_sigma2") else {
    c("samples_U","samples_A","samples_tau","samples_logtheta","samples_sigma2")
  }
  max(vapply(names,function(name) {
    x<-a$mcmc[[name]];y<-b$mcmc[[name]]
    if(!identical(dim(x),dim(y)) || length(x)!=length(y) ||
       any(!is.finite(x)) || any(!is.finite(y))) return(Inf)
    max(abs(x-y))
  },numeric(1)))
}

summary_rows <- list();comparison_rows <- list();case_id <- 0L
for(study in cfg$studies) for(kernel in cfg$kernels) for(n_calib in cfg$calibration) {
  case_id <- case_id+1L
  e <- engines[[study]]
  data <- datasets[[study]]
  idx <- head(calibration_order,n_calib)
  fit_seed <- cfg$seed+1000L*as.integer(study)+10L*n_calib
  paired <- list()
  ## Alternate which implementation runs first to reduce a systematic warmup bias.
  modes <- if(case_id %% 2L) c(FALSE,TRUE) else c(TRUE,FALSE)
  for(schur in modes) {
    mode <- if(schur) "cache_schur" else "cache_dense"
    key <- paste0("study",study,"_",kernel,"_calib",n_calib,"_",mode)
    cat(format(Sys.time(),"%H:%M:%S")," ",key," ...\n",sep="");flush.console()
    warnings <- character(0)
    started <- proc.time()[[3L]]
    fit <- tryCatch(withCallingHandlers({
      if(study=="1") {
        observed <- rep(NA_real_,cfg$n);observed[idx] <- data$train$u[idx]
        e$fit_eivgp_1d(data$train$x,data$train$y,data$train$c,
          u_obs=observed,calib_idx=idx,m=data$m,n_iter=cfg$n_iter,burn=cfg$burn,
          thin=1L,n_chains=cfg$n_chains,preset=cfg$preset,seed=fit_seed,
          parallel_chains=cfg$parallel,n_cores=cfg$n_cores,kernel=kernel,
          gp_block_schur=schur)
      } else {
        e$fit_eivgp_ordprobit_fb(data$train$X,data$train$y,data$train$C,
          U_obs=data$train$U,calib_idx=idx,d=2L,m_vec=rep(4L,4L),
          n_iter=cfg$n_iter,burn=cfg$burn,thin=1L,n_chains=cfg$n_chains,
          preset=cfg$preset,seed=fit_seed,parallel_chains=cfg$parallel,
          n_cores=cfg$n_cores,kernel=kernel,store_scores=FALSE,
          control_overrides=list(gp_use_block_schur=schur))
      }
    },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=identity)
    elapsed <- proc.time()[[3L]]-started
    if(inherits(fit,"error")) {
      saveRDS(list(error=conditionMessage(fit),warnings=warnings,elapsed=elapsed,
        config=cfg),file.path(out_dir,paste0(key,".rds")))
      stop("Computational fit failed for ",key,": ",conditionMessage(fit),
           ". Completed output remains in ",out_dir)
    }
    sampler_seconds <- fit$diagnostics$summary$time_seconds
    diag <- diagnose(fit,study,sampler_seconds)
    counts <- telemetry(fit)
    saveRDS(list(fit=fit,diagnostics=diag,telemetry=counts,warnings=warnings,
      elapsed_seconds=elapsed,config=cfg,benchmark_only=TRUE),
      file.path(out_dir,paste0(key,".rds")))
    write.csv(diag$detail,file.path(out_dir,paste0(key,"_diagnostics.csv")),row.names=FALSE)
    write.csv(counts,file.path(out_dir,paste0(key,"_telemetry.csv")),row.names=FALSE)
    stat <- diag$selected_stats
    row <- data.frame(study=study,kernel=kernel,n=cfg$n,n_calib=n_calib,
      mode=mode,n_iter=cfg$n_iter,burn=cfg$burn,n_chains=cfg$n_chains,
      sampler_seconds=sampler_seconds,elapsed_seconds=elapsed,
      all_saved_values_finite=finite_fit(fit,study),
      raw_diagnostics_finite=diag$raw_stats$finite,
      raw_max_rhat=diag$raw_stats$max_rhat,
      raw_min_bulk_ess=diag$raw_stats$min_bulk,raw_min_tail_ess=diag$raw_stats$min_tail,
      raw_strict_pass=diag$raw_pass,selected_scope=diag$scope,
      selected_diagnostics_finite=stat$finite,selected_max_rhat=stat$max_rhat,
      selected_min_bulk_ess=stat$min_bulk,selected_min_tail_ess=stat$min_tail,
      selected_min_bulk_ess_per_second=stat$min_bulk/sampler_seconds,
      selected_min_tail_ess_per_second=stat$min_tail/sampler_seconds,
      selected_strict_pass=diag$selected_pass,
      warning_count=length(warnings),publication_evidence=FALSE,
      stringsAsFactors=FALSE)
    summary_rows[[length(summary_rows)+1L]] <- row
    paired[[mode]] <- list(fit=fit,row=row,diagnostics=diag)
    write.csv(do.call(rbind,summary_rows),file.path(out_dir,"fit-summary.csv"),row.names=FALSE)
    cat("  ",round(sampler_seconds,3)," s; selected Rhat ",round(stat$max_rhat,3),
        "; min bulk ESS ",round(stat$min_bulk,1),"; strict diagnostic pass ",
        diag$selected_pass,"\n",sep="");flush.console()
  }
  dense <- paired$cache_dense;schur <- paired$cache_schur
  difference <- path_difference(dense$fit,schur$fit,study)
  comparison_rows[[length(comparison_rows)+1L]] <- data.frame(
    study=study,kernel=kernel,n_calib=n_calib,
    dense_seconds=dense$row$sampler_seconds,schur_seconds=schur$row$sampler_seconds,
    dense_over_schur_speed_ratio=dense$row$sampler_seconds/schur$row$sampler_seconds,
    max_saved_path_difference=difference,
    same_random_path_within_tolerance=is.finite(difference)&&difference<=cfg$path_tolerance,
    selected_bulk_ess_per_second_ratio=schur$row$selected_min_bulk_ess_per_second/
      dense$row$selected_min_bulk_ess_per_second,
    selected_tail_ess_per_second_ratio=schur$row$selected_min_tail_ess_per_second/
      dense$row$selected_min_tail_ess_per_second,
    both_selected_strict_pass=dense$row$selected_strict_pass&&schur$row$selected_strict_pass,
    publication_evidence=FALSE,stringsAsFactors=FALSE)
  write.csv(do.call(rbind,comparison_rows),file.path(out_dir,"paired-comparison.csv"),row.names=FALSE)
}
cat("\nBenchmark comparisons (diagnostic failure is reported, never hidden):\n")
print(do.call(rbind,comparison_rows),row.names=FALSE)
cat("\nSaved fits, strict diagnostics, source hashes and counters in ",out_dir,"\n",sep="")
cat("This benchmark does not evaluate every scientific functional and is not a publication release gate.\n")
