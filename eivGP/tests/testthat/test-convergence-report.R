test_that("failure attribution and screened matching preserve denominators", {
root<-tempfile();dir.create(root);out<-tempfile()
jobs<-expand.grid(rep=1:3,method=c('UC-GP-CC','LVGP-CC'),stringsAsFactors=FALSE);jobs$calib<-20L;jobs$job<-paste(jobs$method,jobs$rep,sep='-')
saveRDS(jobs,file.path(root,'job_registry.rds'))
for(i in seq_len(nrow(jobs))) {
 j<-jobs[i,];f<-file.path(root,j$job);dir.create(f)
 shared<-file.path(root,sprintf('shared-rep%03d-cal%02d',j$rep,j$calib));dir.create(shared,showWarnings=FALSE)
 write.csv(data.frame(rhat=if(j$rep==2)1.2 else 1,ess_bulk=500,ess_tail=500),file.path(shared,'measurement_diagnostics.csv'),row.names=FALSE)
 if(j$method=='LVGP-CC'&&j$rep==3){saveRDS(list(message='native failure',stage='optimizer'),file.path(f,'failure.rds'));next}
 saveRDS(list(optimization_status=if(j$method=='UC-GP-CC'&&j$rep==3)'package_selected_nonconverged'else'converged'),file.path(f,'predictions.rds'))
 saveRDS(data.frame(target='y',stratum='seen',coordinate=NA_integer_,sq_error=if(j$method=='UC-GP-CC')1 else 4,crps=1,nlpd=1,coverage=1,width=2),file.path(f,'scores.rds'))
}
z<-report_numerical_convergence(c(test=root),out)
stopifnot(sum(z$status$converged)==2,all(z$matched_converged_summary$n_matched==1),all(z$failure_rates$measurement_n_failed==1),all(z$failure_rates$optimizer_n_failed==1))
stopifnot(all(z$matched_converged_summary$difference[z$matched_converged_summary$metric=='rmse']==1)) # LVGP minus UC-GP
unlink(file.path(root,jobs$job[1],'predictions.rds'))
z<-report_numerical_convergence(c(test=root),out);stopifnot(!nrow(z$matched_converged_summary),z$status$optimizer[1]=='unknown')
stopifnot(inherits(try(report_numerical_convergence(c(test=root),file.path(root,'report')),silent=TRUE),'try-error'))
stopifnot(inherits(try(report_numerical_convergence(c(test=root,test=root),out),silent=TRUE),'try-error'))
writeLines("corrupt",file.path(root,jobs$job[1],'diagnostics.rds'))
z<-report_numerical_convergence(c(test=root),out)
expect_equal(nrow(z$read_errors),1L)
expect_identical(z$status$optimizer[1],'unknown')
unlink(c(root,out),recursive=TRUE)
cat('Failure attribution, unknown exclusion, denominators, matched pairing, duplicate rejection and source protection passed.\n')

})
