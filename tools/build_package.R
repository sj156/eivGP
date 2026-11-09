# Run from the repository root. Refresh source provenance without rewriting code.
root <- normalizePath('eivGP',mustWork=TRUE)
modules <- c(file.path(root,'R',c('helpers.R','mcmc.R','hmc.R','RcppExports.R','plotting.R','disk.R')),
             list.files(file.path(root,'src'),'([.]cpp$|^Makevars)',full.names=TRUE))
definitions <- character()
for(path in list.files(file.path(root,'R'),'[.]R$',full.names=TRUE)) {
  expressions <- parse(path)
  for(expression in expressions) {
    if(!is.call(expression) || !is.symbol(expression[[1L]]) ||
       !as.character(expression[[1L]])%in%c('<-','=') ||
       !is.symbol(expression[[2L]]) || !is.call(expression[[3L]]) ||
       !identical(expression[[3L]][[1L]],as.name('function')))next
    name <- as.character(expression[[2L]])
    if(name%in%definitions)stop('Duplicate package function: ',name,' in ',path)
    definitions <- c(definitions,name)
  }
}
utils::write.csv(data.frame(module=basename(modules),md5=unname(tools::md5sum(modules))),
                file.path(root,'inst','module-manifest.csv'),row.names=FALSE)
cat('R and native sampler sources parsed/hashed. Install with R CMD INSTALL eivGP.\n')
