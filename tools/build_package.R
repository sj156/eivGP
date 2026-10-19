# Run from the repository root. Refresh source provenance without rewriting code.
root <- normalizePath('eivGP',mustWork=TRUE)
modules <- c(file.path(root,'R',c('helpers.R','mcmc.R','hmc.R','RcppExports.R','plotting.R')),
             list.files(file.path(root,'src'),'([.]cpp$|^Makevars)',full.names=TRUE))
for(path in list.files(file.path(root,'R'),'[.]R$',full.names=TRUE))parse(path)
utils::write.csv(data.frame(module=basename(modules),md5=unname(tools::md5sum(modules))),
                file.path(root,'inst','module-manifest.csv'),row.names=FALSE)
cat('R and native sampler sources parsed/hashed. Install with R CMD INSTALL eivGP.\n')
