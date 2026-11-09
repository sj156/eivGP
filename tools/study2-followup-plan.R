# Study II targeted continuation plan. Sourcing or executing this file never fits.
# Sampling uses eivGP::continue_eivgp_disk(); full-history assessment is a separate phase.
# tools/study2.R reports saved results; the package Study II driver calls the ESS
# module when STUDY2_ESS_MIN is set. The fixed-budget paper lane sets it to NULL.

study2_followup_plan <- function(repo_root, source_run=NULL, output_root=NULL,
                                cores=12L, scope=c('pilot','all'),data_root=NULL) {
  scope<-match.arg(scope)
  repo_root<-normalizePath(repo_root,mustWork=TRUE)
  if(is.null(source_run)) {
    candidates<-list.dirs(file.path(repo_root,'results/study2-eivGP/core'),recursive=FALSE,full.names=TRUE)
    candidates<-candidates[file.exists(file.path(candidates,'combined/mcmc_diagnostics.csv'))]
    if(length(candidates)!=1L)stop('Supply --source=PATH to one saved Study II run with combined/mcmc_diagnostics.csv')
    source_run<-candidates[1L]
  }
  if(is.null(output_root))output_root<-file.path(repo_root,'results-sensitivity/study2-shared-continuation-v1')
  if(is.null(data_root)) {
    local_data<-file.path(repo_root,'results/data/study2')
    data_root<-if(dir.exists(local_data))local_data else file.path(path.expand('~/eivGP-results'),'data/study2')
  }
  canonical<-function(p) {
    p<-path.expand(p)
    if(!grepl('^/',p))p<-file.path(repo_root,p)
    if(file.exists(p))return(normalizePath(p,mustWork=TRUE))
    file.path(canonical(dirname(p)),basename(p))
  }
  source_run<-canonical(source_run);output_root<-canonical(output_root);data_root<-canonical(data_root)
  protected<-c(canonical(file.path(repo_root,'results')),source_run,data_root)
  if(any(vapply(protected,function(p)output_root==p||startsWith(output_root,paste0(p,'/'))||startsWith(p,paste0(output_root,'/')),logical(1))))
    stop('The follow-up output must be separate from the existing results tree')
  if(length(cores)!=1L||!is.numeric(cores)||!is.finite(cores)||cores<4L||cores!=as.integer(cores))stop('cores must be an integer >=4')
  d<-read.csv(file.path(source_run,'combined/mcmc_diagnostics.csv'),stringsAsFactors=FALSE)
  d<-d[d$cell_id=='primary_q4_calibration'&d$n_calib%in%c(0L,20L),,drop=FALSE]
  if(nrow(d)!=100L||anyDuplicated(d[c('rep','n_calib')]))stop('Expected 100 unique q=4 calibration-0/20 source fits')
  pilot<-c(5L,8L,16L,25L,42L)
  if(scope=='pilot')d<-d[d$rep%in%pilot,,drop=FALSE]
  d<-d[order(d$rep,d$n_calib),,drop=FALSE]
  jobs<-data.frame(rep=as.integer(d$rep),calibration=as.integer(d$n_calib),
    original_gate_pass=d$mcmc_pass,scientific_pass=d$target_functional_pass,
    invariant_pass=d$invariant_measurement_pass,raw_pass=d$raw_coordinate_pass)
  jobs$fit<-mapply(function(rep,cal) {
    p<-list.files(file.path(source_run,'cells/primary_q4_calibration'),
      sprintf('^recovery_fit_primary_rep%03d_cal%03d_.*[.]rds$',rep,cal),recursive=TRUE,full.names=TRUE)
    if(length(p)!=1L)stop('Missing or ambiguous original fit: ',rep,'/',cal)
    normalizePath(p)
  },jobs$rep,jobs$calibration)
  jobs$data<-file.path(data_root,'primary_q4_calibration',
    sprintf('study2_primary_q4_n100_ntest100_m4_rep%03d.rds',jobs$rep))
  if(any(!file.exists(jobs$data)))stop('Missing frozen datasets')
  jobs$output<-file.path(output_root,sprintf('rep%03d-cal%02d',jobs$rep,jobs$calibration))
  structure(list(
    schema=1L,scope=scope,source_run=source_run,output_root=output_root,jobs=jobs,
    resources=list(cores=as.integer(cores),chains_per_fit=4L,
      concurrent_fits=min(nrow(jobs),as.integer(cores)%/%4L),blas_threads=1L),
    implementation=list(
      continuation_api='eivGP::continue_eivgp_disk',
      module='eivGP/R/disk.R',driver='eivGP/R/helpers.R',
      report_entrypoint='tools/study2.R',
      fitting_entrypoint='tools/study2-eivGP.R',
      continuation_entrypoint='tools/study2-continue-eivGP.R',
      reuse='Use the existing disk continuation module; do not load or concatenate retained history in sampling workers.',
      paper_lane='Keep fixed-budget results as published; do not enable continuation in their output root.'),
    acceptance=list(rhat_max=1.01,bulk_ess_min=400L,tail_ess_min=400L,
      mcse_over_posterior_sd_max=.05,min_chains=4L,full_retained_history=TRUE,
      required=c('scientific_targets','ordinal_measurement_invariants',
        'raw_free_coordinates_when_calibration_positive'),
      undefined_required_diagnostic='fail',trace_review_required=TRUE),
    phases=c('serial_legacy_import','bounded_disk_sampling','serial_full_history_assessment'),
    sampling=list(chunk_iter=500L,default_additional_iterations=2000L,
      diagnostics_during_sampling=FALSE,automatic_ess_extensions=FALSE),
    shared_module_call=quote(eivGP::continue_eivgp_disk(
      path=archive_path,n_iter=stage_iterations,chunk_iter=500L,
      parallel=TRUE,n_cores=4L)),
    module_contract=c(
      'Copy/convert each source fit into a separate output root before importing it with n_iter=0.',
      'Import is serial and reads one legacy fit once; sampling uses only terminal states, indexes and bounded buffers.',
      'All retained draws and HMC history remain in immutable disk blocks; keep fit.rds with fit.rds.draws.',
      'Each explicit request has a fixed target so repeating it resumes without extending twice.',
      'Run full-history diagnostics separately after sampling. read_fit() materializes history and must never be called by the MCMC executor.'),
    contingencies=c(
      baseline='Validate fit/data/calibration identities, assess existing draws, freeze scientific panel and seeds, and save source hashes.',
      low_ess_stable_chains='After a separate assessment, request a short explicit disk-backed extension; no diagnostic work runs during MCMC.',
      persistent_disagreement='If offline diagnostics show separated chains, do not issue further sampling requests before identification/initialization/sampler review.',
      adequate_ess_other_failure='Scientific ESS success alone is insufficient; inspect invariants/raw coordinates, target uncertainty and MCSE before another bounded extension.',
      passed_pilot='Review complete traces and scientific-panel coverage, then use scope=all to include all 100 datasets; retain failures in comparisons.',
      exhausted='Stop at the requested bounded target. Assess after sampling before explicitly requesting another round.',
      separated_modes='Consider independent four-chain refits with dispersed valid initial states on the same frozen data, priors and likelihood, in another labeled root. Do not pool unresolved modes or discard unfavorable history.',
      scoring='After diagnostic review, recompute original and continued scores using identical settings, check integration-seed/draw sensitivity, and compare paired changes against replication MCSE.'),
    preservation=c(
      'Never update original fits, reports, caches or checkpoints; never run the whole paper workflow in the original root to continue these fits.',
      'Use an isolated installed package and frozen source manifest. Do not replace a library used by any active job.',
      'The earlier long pilot is stopped. Its immutable saved states can be copied into a new disk archive; never migrate the originals in place.',
      'This plan does not guarantee every check will pass; unresolved cases remain visible.'),
    execution_enabled=FALSE),class='study2_followup_plan')
}

if(sys.nframe()==0L) {
  if('--help'%in%commandArgs(TRUE)) {
    cat('Plan only: --source=RUN_ROOT --data-root=DATA_ROOT --output=NEW_ROOT --scope=pilot|all --cores=12\n')
    quit(status=0L)
  }
  opts<-list(repo=NULL,source=NULL,output=NULL,scope='pilot',cores='12',`data-root`=NULL)
  for(arg in commandArgs(TRUE)) {
    if(!grepl('^--[^=]+=.+$',arg))stop('Use --name=value')
    key<-sub('^--([^=]+)=.*$','\\1',arg)
    if(!key%in%names(opts))stop('Unknown option: ',key)
    opts[[key]]<-sub('^--[^=]+=','',arg)
  }
  if(is.null(opts$repo)) {
    script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1L])
    opts$repo<-dirname(dirname(normalizePath(script)))
  }
  plan<-study2_followup_plan(opts$repo,opts$source,opts$output,as.numeric(opts$cores),opts$scope,opts[['data-root']])
  print(plan$jobs[c('rep','calibration','original_gate_pass')],row.names=FALSE)
  print(plan$resources);print(plan$phases);print(plan$sampling);print(plan$shared_module_call)
  cat('\nPlan only: no fitting, continuation, installation, or existing-result updates.\n')
  # Explicit --output requests a separate, reviewable plan artifact.
  if(!is.null(opts$output)) {
    dir.create(plan$output_root,recursive=TRUE,showWarnings=FALSE)
    dest<-file.path(plan$output_root,paste0('followup-plan-',plan$scope,'.rds'))
    manifest<-file.path(plan$output_root,paste0('plan-',plan$scope,'.csv'))
    if(any(file.exists(c(dest,manifest))))stop('Plan artifact already exists; choose a fresh output directory')
    saveRDS(plan,dest);write.csv(plan$jobs,manifest,row.names=FALSE)
    cat('Saved: ',dest,'\nManifest: ',manifest,'\n')
  }
}
