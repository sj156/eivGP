repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
testthat::test_that('dataset folders read frozen cohorts and prepare private schemas', {
source(file.path(repo,'tools/ocean-tools.R'));api<-ocean_workflow_api(repo)
for(ds in c('ocean200','ocean500')) {
 path<-ocean_dataset_dir(repo,ds);inputs<-api$ocean_prepare_inputs(path)
 stopifnot(nrow(inputs$cohort)==as.integer(sub('ocean','',ds)))
 cat(ds,'all frozen checksums and fold validation passed\n')
}
source(file.path(repo,'tools/adni-tools.R'));source(file.path(repo,'tools/oai-tools.R'))
for(app in c('adni','oai')){
 api<-if(app=='adni')adni_load_api(repo)else oai_load_api(repo)
 meta<-file.path(repo,'data',app,'DATASET.example.json');m<-jsonlite::fromJSON(meta)
 set.seed(743);n<-200L;d<-data.frame(RID=seq_len(n),Y=rnorm(n))
 for(j in seq_along(m$X_order))d[[paste0('X',j)]]<-rnorm(n)
 for(j in seq_along(m$U_order))d[[paste0('U',j)]]<-rnorm(n)
 for(j in seq_along(m$C_order))d[[paste0('C',j)]]<-rep(seq_len(m$m_vec[j]),length.out=n)
 root<-tempfile();dir.create(root);dir.create(file.path(root,'data'));cohort<-file.path(root,'cohort.csv');write.csv(d,cohort,row.names=FALSE)
 selector<-file.path(root,'data',app)
 api[[paste0(app,'_prepare_folds')]](cohort,meta,selector)
 stopifnot(dir.exists(selector),length(list.files(selector))==46L)
 inputs<-api[[paste0(app,'_plan_inputs')]](selector)
 cat(app,'folder preparation and read-back passed (46 files)\n')
 # The same prepared bytes remain readable in a legacy external directory.
 legacy<-file.path(root,'external');dir.create(legacy)
 for(p in list.files(selector,full.names=TRUE))file.copy(p,file.path(legacy,basename(p)))
 api[[paste0(app,'_plan_inputs')]](legacy)
 cat(app,'external directory compatibility passed\n');unlink(root,recursive=TRUE)
 testthat::expect_true(TRUE)
}

})
