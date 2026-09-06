# setwd('C:/Users/pengh/OneDrive/Research/Projects/Mixed-input GP/Experiments')
setwd('C:/Users/penghui/OneDrive/Research/Projects/Mixed-input GP/Experiments')

library(mvtnorm)
library(TruncatedNormal)
library(scatterplot3d)
library(lhs)
library(calculus)

source('Basics.R')

# Generating Data -----
set.seed(6)
category.type <- "ordinal"  
m <- 3  # number of levels of c
nc <- 1  # dimension of x
D <- ifelse(category.type == 'ordinal',1, m-1) # dimension of u

if (category.type == 'ordinal'){  # generating Tau
  Tau <- sort(rnorm(m-1))  
}else{
  Tau <- rnorm(m-1)  
}

n <- 100  # Data size
XU <- 4*maximinLHS(n, nc + D) - 2 # generating x and u by maximin LHS on (-2,2)
X <- as.matrix(XU[, 1:nc]) 
U <- as.matrix(XU[, (nc+1):(nc+D)])
W <- cbind(X, U)
c <- apply(U, 1, function(u) match(T, c(u-Tau,-Inf) <= 0))  # generating c

coef1 <- 1
coef2 <- 1
mean.term <- MeanFunction(X, U, coef1=coef1, coef2=coef2)  # the mean term

sigma <- 0.1 # the magnitude of the Gaussian process
sigma.error <- 0.1 # the std of error term
length.scale.x <- 1
length.scale.u <- 1
Theta <- c(rep(length.scale.x, nc), rep(length.scale.u, D), sigma^2) # GP hyperparameters, including length scales and variance
Kn <- KernelMatrix(W, Theta) 
Cn <- Kn + sigma.error^2 * diag(n) 
y <- mean.term + c(rmvnorm(1, rep(0, n), Cn))  
#y <- mean.term + rnorm(n,0,sigma.error)  
# generating y

# Sampling from the posterior by Gibbs ----
initial.values <- GenerateInitial(3, Tau0 = c(0, 1),
                                  category.type = category.type, n = n, D = D, m = m, c = c)
# creating empty matrices to store samples
U.hist <- initial.values$U0
if (category.type == 'ordinal')
  U.hist <- matrix(U.hist, ncol = 1) # avoid dropping
Tau.hist <- matrix(initial.values$Tau0, nrow = m-1)
accep.hist <- NULL

K = 1000
time.hist <- c(timestamp())
lmd <- 0.05
stepsize.hist <- data.frame(stepsize = lmd, k = 1)
for (cycles in 1:60) {
  Gibbs1 <- Gibbs(seed = NULL, Tau0 = Tau.hist[, ncol(Tau.hist)], U0 = U.hist[, (ncol(U.hist)-D+1):ncol(U.hist)], K = K, Metropolis = 'Langevin', lmd = lmd,
                  n = n, D = D, m = m, 
                  c = c, X = X, y = y, 
                  Theta = Theta, sigma.error = sigma.error,
                  consider.mean = T, coef1 = as.matrix(coef1), coef2 = as.matrix(coef2))
  U.hist <- cbind(U.hist, Gibbs1$U.hist)
  Tau.hist <- cbind(Tau.hist, Gibbs1$Tau.hist)
  accep.hist <- c(accep.hist, Gibbs1$accep)
  if (mean(Gibbs1$accep) <= 0.2) {
    lmd <- lmd / 1.5
    stepsize.hist[nrow(stepsize.hist) + 1, ] <- c(lmd, K*cycles+1)
  }else if (mean(Gibbs1$accep) > 0.5){
    lmd <- lmd * 1.5
    stepsize.hist[nrow(stepsize.hist) + 1, ] <- c(lmd, K*cycles+1)
  }
  time.hist <- c(time.hist, timestamp())
  save.image(file = 'MCMC.results.RData')
}

system('shutdown -s')



