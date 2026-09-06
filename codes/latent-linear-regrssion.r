
datobj <- data.frame(y=1:4, c=c(rep(0,2),rep(1,2)))

plot(datobj$c, datobj$y)


n <- length(datobj$y)
Nmc <- 1e3
p <- 2 ## linear regression with intercept 

### prior
nu0 <- 3
S0 <- 1 
### store posterior samples 
pos.u <- matrix(NULL, nrow = n, ncol = Nmc)
pos.beta <- matrix(NULL, nrow = p, ncol = Nmc)
pos.sig <- matrix(NULL, nrow = 1, ncol = Nmc)

pos.tau <- matrix(NULL, nrow = 1, ncol = Nmc)
pos.sigu <- matrix(NULL, nrow = 1, ncol = Nmc) 

## initialize 
beta <- rnorm(p)
sig <- 1
tau <- 0  

for(ii in 1:Nmc){
    uu <- 
}
