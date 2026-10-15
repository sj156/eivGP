test_that("identity conflicts explain sampling differences and preserve saved results", {
  output<-tempfile("ocean-saved-");dir.create(output)
  on.exit(unlink(output,recursive=TRUE),add=TRUE)
  saved<-list(data=c(input="abc"),code="old",packages=c(eivGP="0.4.0.9007"),
    numerics=list(chains=4L,iter=100000L,burn=5000L))
  current<-saved;current$numerics$iter<-50000L;current$code<-"new"
  saveRDS(saved,file.path(output,"workflow_identity.rds"))
  saveRDS(list(draws=1:10),file.path(output,"fit.rds"))
  dir.create(paste0(output,"_new"))
  on.exit(unlink(paste0(output,"_new"),recursive=TRUE),add=TRUE)
  files<-list.files(output,full.names=TRUE)
  hashes<-tools::md5sum(files);times<-file.info(files)$mtime
  message<-tryCatch(ocean_assert_identity(saved,current,output),error=conditionMessage)
  expect_match(message,"Changed identity fields: code, numerics",fixed=TRUE)
  expect_match(message,"iter: saved=100000; requested=50000",fixed=TRUE)
  expect_match(message,paste0("--output=",shQuote(paste0(output,"_new2"))),fixed=TRUE)
  expect_match(message,"Do not delete workflow_identity.rds",fixed=TRUE)
  expect_false(dir.exists(paste0(output,"_new2")))
  expect_identical(tools::md5sum(files),hashes)
  expect_identical(file.info(files)$mtime,times)
  expect_true(ocean_assert_identity(saved,saved,output))
  changed<-saved;changed$packages<-c(eivGP="0.4.0.9008")
  expect_error(ocean_assert_identity(saved,changed,output),"Changed identity fields: packages")
  changed<-saved;changed$data<-c(input="changed")
  expect_error(ocean_assert_identity(saved,changed,output),"Changed identity fields: data")
})

test_that("the executor stops incompatible runs before any task starts", {
  output<-tempfile("ocean-conflict-");dir.create(output)
  on.exit(unlink(output,recursive=TRUE),add=TRUE)
  saved<-list(numerics=list(chains=4L,iter=100000L,burn=5000L))
  current<-saved;current$numerics$iter<-50000L
  saveRDS(saved,file.path(output,"workflow_identity.rds"))
  saveRDS(list(draws=1:10),file.path(output,"fit.rds"))
  files<-list.files(output,full.names=TRUE);hashes<-tools::md5sum(files)
  testthat::local_mocked_bindings(
    ocean_prepare_inputs=function(...)list(),
    ocean_select_tasks=function(...)data.frame(module="eiv"),
    ocean_identity=function(...)current,
    ocean_schedule=function(...)stop("Fitting must not start"),.package="eivGP")
  expect_error(ocean_execute("unused","unused",output,ocean_settings(action="fit")),
    "iter: saved=100000; requested=50000",fixed=TRUE)
  expect_identical(tools::md5sum(files),hashes)
  expect_identical(sort(list.files(output,all.files=TRUE,no..=TRUE)),sort(basename(files)))
})

test_that("execution core changes do not invalidate saved fit identities", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE),add=TRUE)
  a<-ocean_identity(root,root,ocean_settings(cores=4))
  b<-ocean_identity(root,root,ocean_settings(cores=16))
  expect_true(ocean_assert_identity(a,b,root))
})


test_that("initial CLI fitting selects a new folder without changing old checkpoints", {
  root<-tempfile("ocean-reroute-");dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  old<-file.path(root,"paper");dir.create(old)
  id<-list(numerics=list(iter=100000L));saveRDS(id,file.path(old,"workflow_identity.rds"))
  saveRDS(1:10,file.path(old,"fit.rds"));files<-list.files(old,full.names=TRUE)
  hashes<-tools::md5sum(files);times<-file.info(files)$mtime
  repo<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  env<-new.env(parent=globalenv());sys.source(file.path(repo,"paper/applications/ocean.R"),env)
  expect_identical(env$ocean_choose_fit_output(old,id),old)
  message<-capture.output(new<-env$ocean_choose_fit_output(old,list(numerics=list(iter=50000L))))
  expect_true(dir.exists(new));expect_false(identical(new,old))
  expect_identical(list.files(new),character())
  expect_true(any(grepl("Existing results preserved",message,fixed=TRUE)))
  expect_true(any(grepl(paste0("--output=",shQuote(new)),message,fixed=TRUE)))
  expect_identical(tools::md5sum(files),hashes);expect_identical(file.info(files)$mtime,times)
  env$ocean_default_output<-function(...)old
  env$ocean_resolve_path<-getFromNamespace("ocean_resolve_path","eivGP")
  saveRDS(new,paste0(old,".latest-run.rds"))
  expect_identical(env$ocean_output_dir(repo,"ocean500","paper"),new)
  expect_identical(env$ocean_output_dir(repo,"ocean500","paper",old),old)
  empty<-file.path(root,"empty")
  expect_identical(env$ocean_choose_fit_output(empty,id),empty)
  expect_false(dir.exists(empty))
  legacy<-file.path(root,"legacy");dir.create(legacy);saveRDS(1,file.path(legacy,"old.rds"))
  invisible(capture.output(legacy_new<-env$ocean_choose_fit_output(legacy,id)))
  expect_true(dir.exists(legacy_new));expect_true(file.exists(file.path(legacy,"old.rds")))
})
