test_that("replication and chain allocation respects the total budget", {
  for(n in c(1L,2L,3L,4L,5L,8L,56L)) {
    a <- mixedgp_followup_allocation(n,50L)
    expect_lte(a$workers*a$chains,n)
    expect_lte(a$chains,4L)
  }
  if(.Platform$OS.type=="unix")expect_equal(mixedgp_followup_allocation(56,50)$workers,14L)
  expect_equal(mixedgp_followup_allocation(56,1)$workers,1L)
})

test_that("independent workers preserve explicit random streams and cached outputs", {
  root <- tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  groups <- as.list(1:6)
  worker <- function(i) {
    p <- file.path(root,paste0(i,".rds"))
    if(!file.exists(p)) {set.seed(100+i);saveRDS(runif(4),p)}
    sum(readRDS(p))
  }
  first <- mixedgp_followup_apply(groups,worker,mixedgp_followup_allocation(1,6))
  files <- list.files(root,full.names=TRUE); before <- tools::md5sum(files)
  second <- mixedgp_followup_apply(groups,worker,mixedgp_followup_allocation(8,6))
  expect_identical(first,second)
  expect_identical(before,tools::md5sum(files))
  unlink(files)
  third <- mixedgp_followup_apply(groups,worker,mixedgp_followup_allocation(8,6))
  expect_identical(first,third)
})

test_that("resume compatibility permits only the recorded scheduler edit", {
  compat <- data.frame(field=c("helpers.R","mechanism"),old=c("oldh","oldm"),new=c("newh","newm"))
  old <- list(package_manifest=c(helpers.R="oldh",mcmc.R="unchanged"),code=c(mechanism="oldm"),settings=list(iter=100))
  new <- old;new$package_manifest["helpers.R"] <- "newh";new$code["mechanism"] <- "newm"
  expect_true(mixedgp_followup_identity_equal(old,new,compat))
  new$settings$iter <- 101
  expect_false(mixedgp_followup_identity_equal(old,new,compat))
  new$settings$iter <- 100;new$package_manifest["mcmc.R"] <- "changed"
  expect_false(mixedgp_followup_identity_equal(old,new,compat))
})
