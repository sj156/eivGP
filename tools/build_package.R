# Run from the repository root. Refresh provenance without rewriting source.
root <- normalizePath("eivGP",mustWork=TRUE)
modules <- file.path(root,"R",c("helpers.R","mcmc.R","plotting.R"))
for(path in list.files(file.path(root,"R"),"[.]R$",full.names=TRUE)) parse(path)
utils::write.csv(data.frame(module=basename(modules),md5=unname(tools::md5sum(modules))),
  file.path(root,"inst","module-manifest.csv"),row.names=FALSE)
cat("Package sources parsed and module manifest updated. Install with R CMD INSTALL eivGP.\n")
