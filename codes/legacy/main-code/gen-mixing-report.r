

dir.output <- "~/Dropbox/code/output/MCMC-image/"
dir.report <- "~/Dropbox/code/output/mixing-report/"
dir.code <- "~/Dropbox/code/real-data-code/" 

setwd(dir.code)

## configs <- cbind(10, paste0('2410',13:18))
## configs <- rbind(configs,
##                  cbind(11, paste0('2410',22:24)))

configs <- cbind(10, paste0('24-11-23-3less')) 

for(ii in 1:nrow(configs)){
    params <- list(numclust =configs[ii,1], timestamp = configs[ii,2])
    output.filename <- paste0('MCMC-mixing-report',
                              params$numclust,'-',
                              params$timestamp,'.pdf')
    
    rmarkdown::render('BMoE-check-mixing.Rmd', envir = new.env() ,
                      params= params,
                      output_file = paste0(dir.report,output.filename))
    
}


