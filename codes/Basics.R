######################################
#############  functions  ############
######################################

# Requiring package: mvtnorm, TruncatedNormal

MeanFunction <- function (X, U, coef1, coef2) {
  # This function calculates a mean function given data matrices X and U.
  return(0.2*sin(pi*(X%*%coef1)) +U%*%coef2)
}

PlotMeanFunction <- function (nx, coef1, coef2) {
  # Draw a perspective plot for the mean function. Only univariate x and u are allowed. 
  w <- seq(-1, 1, length=nx)
  W <- expand.grid(w, w)
  y <- MeanFunction(matrix(W[,1],ncol=1), matrix(W[,2],ncol=1), coef1, coef2)
  persp(w, w, matrix(y, ncol=nx), theta=-30, phi=30, xlab="x",
        ylab="u", zlab="y")
}

LogDensityRatio <- function(y, mu1, sigma1, mu2, sigma2, sigma1i = NULL, sigma2i = NULL) {
  # Calculate the log density ratio of two mtv normal at the same location y
  if(is.null(sigma1i))
    sigma1i <- solve(sigma1)
  if(is.null(sigma2i))
    sigma2i <- solve(sigma2)
  logdet1 <- c(determinant(sigma1)$modulus)
  logdet2 <- c(determinant(sigma2)$modulus)
  if (is.na(logdet1 - logdet2)) {
    stop('Numerical issues when calculating the difference of log determinants of covariance matrices!')
  }
  return(-(t(y-mu1)%*%sigma1i%*%(y-mu1) - t(y-mu2)%*%sigma2i%*%(y-mu2) + logdet1 - logdet2)/2)
}

KernelMatrix <- function(W, Theta) {
  # Computes the kernel matrix of W with the inverse exponential squared kernel.
  #
  # Args:
  #   W: The data matrix with each row being an observation.
  #   Theta: The last one is the overall variance, and the rest are length scales.
  #
  # Returns:
  #   A kernel matrix.
  
  if (length(Theta) != ncol(W)+1) {
    stop("Wrong number of parameters!")
  }
  if (any(Theta <= 0)) {
    stop('Negative scale parameters!')
  }
  
  n <- nrow(W)
  K <- matrix(1, n, n)  # We first compute the correlation matrix (i.e., with the overall variance as 1), whose diagonal elements are always 1.
  theta <- Theta[1:ncol(W)]
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      K[i, j] <- exp( - sum(theta * (W[i, ] - W[j, ])^2))
      K[j, i] <- K[i, j]
    }
  }
  
  K <- Theta[ncol(W)+1] * K  # Finally, scaling the kernel matrix 
  return(K)
}

GenerateInitial <- function(seed, Tau0 = NULL,
                            category.type, n, D, m, c) {
  # Generates initial values of U and Tau CONSISTENT(!) with the observed c. We first simulate Tau, then use a constrained sampling for U given Tau and c similar to the Gibbs. 
  #
  # Args:
  #   seed: random seed
  #
  # Returns:
  #   A list of initial values of U and Tau.
  
  
  set.seed(seed)
  
  if (category.type == 'ordinal') {  # generating Tau
    if(!is.null(Tau0)){
      Tau0 <- Tau0
    }else{
      Tau0 <- sort(rnorm(m-1)) 
    } 
    Tau.aug <- c(-Inf, Tau0, Inf)
    U0 <- rtmvnorm(1, rep(0,n), diag(n), Tau.aug[c], Tau.aug[c+1])
  }else{
    if(!is.null(Tau0)){
      Tau0 <- Tau0
    }else{
      Tau0 <- rnorm(m-1)  
    }
    Tau.matrix <- matrix(rep(Tau0, n), nrow = n, byrow = T)
    ind.matrix <- matrix(rep(1:D, n), nrow = n, byrow = T)
    ind.matrix.c <- matrix(rep(c, D), ncol = D)
    ind.matrix.lb <- ind.matrix - ind.matrix.c < 0  # logical matrix indicating the indices of U with lower bound constraints
    ind.matrix.ub <- ind.matrix - ind.matrix.c == 0  # logical matrix indicating the indices of U with upper bound constraints (at most 1 True in each row)
    n.lb <- sum(ind.matrix.lb)
    n.ub <- sum(ind.matrix.ub)
    U0 <- matrix(rnorm(n*D), nrow = n)
    U0[ind.matrix.lb] <- rtmvnorm(1, rep(0, n.lb), diag(n.lb), Tau.matrix[ind.matrix.lb], rep(Inf, n.lb))
    U0[ind.matrix.ub] <- rtmvnorm(1, rep(0, n.ub), diag(n.ub), rep(-Inf, n.ub), Tau.matrix[ind.matrix.ub])
  }
  
  return(list(Tau0 = Tau0, U0 = U0))
}

Gradient <- function(U, mu, K, C = NULL, Ci = NULL, length.scale.u = 1) {
  # Calculate the gradient of the log full conditional of latent variable U
  #
  # Args:
  #  U: latent variable 
  #  mu: mean vector at U
  #  K: kernel matrix K_{WW} at U
  #  C: K_{WW} + \sigma_\epsilon^2 I
  #  Ci: C^{-1}. 
  #
  # Return:
  # Gradient (of the same size of U)
  if(is.null(C))
    C <- K + sigma.error^2 * diag(n)
  if(is.null(Ci))
    Ci <- solve(C)
  r <- Ci%*%(y-mu)
  C.tilde <- (Ci - r%*%t(r))*K
  B <- t(matrix(rep(coef2, n), ncol = n))
  D.U <- diag(c(r))%*%B + 2*length.scale.u*(diag(apply(C.tilde, 1, sum)) - C.tilde)%*%U
  return(D.U - U)
}

Hessian <- function(U, mu, K, C = NULL, Ci = NULL, length.scale.u = 1){
  # Calculate the gradient of the log full conditional of latent variable U (ordinal case)
  
  U <- c(U)
  if(is.null(C))
    C <- K + sigma.error^2 * diag(n)
  if(is.null(Ci))
    Ci <- solve(C)
  
  r <- Ci%*%(y-mu)
  U.m.U <- matrix(rep(U, n), ncol = n) - matrix(rep(U, n), ncol = n, byrow = T) # matrix of (ui - uj)
  Q <- K * U.m.U # matrix of Kij(ui-uj) 
  D.r <- Ci%*%(2*length.scale.u*(diag(c(Q%*%r)) + t(Q)%*%diag(c(r))) - c(coef2)*diag(n))
  
  Term.1 <- c(coef2) * t(D.r)
  Term.2 <- 2*length.scale.u*(Ci - r%*%t(r))*(2*length.scale.u*Q*U.m.U - K)
  diag(Term.2) <- 0
  diag(Term.2) <- - apply(Term.2, 1, sum)
  Term.3 <- matrix(0, ncol = n, nrow = n)
  for(s in 1:n) {
    for(l in 1:n) {
      Term.3[s, l] <- 4*length.scale.u^2*(sum(Q[l,]*Ci[,s])*sum(Q[s,]*Ci[,l]) + Ci[s,l]*sum(Q[l,]%*%t(Q[s,])*Ci)) - 2*length.scale.u*(r[l]*sum(Q[l,]*D.r[,s]) + D.r[l,s]*sum(r*Q[l,]))
    }
  }
  
  return(Term.1 + Term.2 + Term.3)
}


Gibbs <- function(seed, Tau0, U0, K, Metropolis, lmd,
                  n, D, m, 
                  c, X, y, 
                  Theta, sigma.error,
                  consider.mean = F, coef1 = NULL, coef2 = NULL) {
  # Runs a a Gibbs sampler for the posterior distribution of Tau and U. 
  #
  # Args:
  #   seed: random seed
  #   Tau0: initial Tau
  #   U0: initial U
  #   K: number of iterations
  #   Metropolis: 'Independent' for the independent prior proposal;
  #               'Random Walk' for the random walk proposal;
  #               'Langevin' for LMC;
  #               'Newton' for Newton MH
  #   lmd: step size (std of the Gaussian proposal)
  #
  # Returns:
  #   A list of MCMC samples of U, Tau, and a binary vector recording acceptances of the Metropolis steps for U. 
  
  if(!is.null(seed))
    set.seed(seed)
  
  U.hist <- matrix(0, nrow = n, ncol = D*K)
  Tau.hist <- matrix(0, nrow = m-1, ncol = K)
  accep <- rep(0, K)
  
  # initializing U, W, KWW, C, Ci, Tau, and mu, (and gradient/Hessian) ----
  Tau <- Tau0  
  U <- U0
  W <- cbind(X, U)
  KWW <- KernelMatrix(W, Theta)
  C <- KWW + sigma.error^2 * diag(n)
  Ci <- solve(C)
  if (consider.mean == F) {
    mu <- rep(0, n)  # zero mean
  }else {
    if (is.null(coef1) || is.null(coef2)) {
      stop('The two coefficients in the mean function are not specified!')
    }
    mu <- MeanFunction(X, U, coef1, coef2)
  }
  if (Metropolis == 'Langevin') {
    GU <- Gradient(U = U, mu = mu, K = KWW, C = C, Ci = Ci)
  }
  if (Metropolis == 'Newton') {
    GU <- Gradient(U = U, mu = mu, K = KWW, C = C, Ci = Ci)
    HU <- Hessian(U = U, mu = mu, K = KWW, C = C, Ci = Ci) - diag(n)
    HUi <- solve(HU)
    U.Newton <- U - lmd^2/2*HUi%*%GU # the proposal mean, i.e., a Newton step
    
    HU1 <<- HU
    HUi1 <<- HUi
    U.Newton1 <<- U.Newton
  }
  
  # Creating some index matrices for subsetting ----
  if (category.type == 'ordinal') {
    ind <- rep(T, m)
    for (j in 1:m){
      ind[j] <- any(c == j)
    }
    ind[1] <- T  
    ind[m] <- T  #For convenience, we can always augment u[c==1] with -Inf and u[c==m] with Inf. Then ind[1] and ind[m] are both T.
  }else{
    ind.matrix <- matrix(rep(1:D, n), nrow = n, byrow = T)
    ind.matrix.c <- matrix(rep(c, D), ncol = D)
    ind.matrix.lb <- ind.matrix - ind.matrix.c < 0  # logical matrix indicating the indices of U with lower bound constraints
    ind.matrix.ub <- ind.matrix - ind.matrix.c == 0  # logical matrix indicating the indices of U with upper bound constraints (at most 1 True in each row)
    n.lb <- sum(ind.matrix.lb) 
    n.ub <- sum(ind.matrix.ub) 
  }
  
  # Gibbs ----
  for(k in 1:K) {
    # Simulating U ----
    # Proposing U.prime ----
    if (category.type == 'ordinal') {
      Tau.aug <- c(-Inf, Tau, Inf)
      if (Metropolis == 'Independent') {
        U.prime <- rtmvnorm(1, rep(0,n), diag(n), Tau.aug[c], Tau.aug[c+1])
      }else if (Metropolis == 'Random Walk') {
        U.prime <- rtmvnorm(1, U, lmd^2 * diag(n), Tau.aug[c], Tau.aug[c+1])
      }else if (Metropolis == 'Langevin') {
        U.prime <- rtmvnorm(1, U + lmd^2/2*GU, lmd^2 * diag(n), Tau.aug[c], Tau.aug[c+1])
      }else if (Metropolis == 'Newton') {
        U.prime <- rtmvnorm(1, U.Newton, -lmd^2 * HUi, Tau.aug[c], Tau.aug[c+1])
      }else{
        stop('Unsupported Metrpolis!')
      }
      
    }else{  # category.type == 'nominal'
      Tau.matrix <- matrix(rep(Tau, n), nrow = n, byrow = T)
      if (Metropolis == 'Independent') {
        U.prime <- matrix(rnorm(n*D), nrow = n)
        U.prime[ind.matrix.lb] <- rtmvnorm(1, rep(0, n.lb), diag(n.lb), Tau.matrix[ind.matrix.lb], rep(Inf, n.lb))
        U.prime[ind.matrix.ub] <- rtmvnorm(1, rep(0, n.ub), diag(n.ub), rep(-Inf, n.ub), Tau.matrix[ind.matrix.ub])
      }else {
        if (Metropolis == 'Random Walk') {
          U.prime.mean <-  U
        }
        if (Metropolis == 'Langevin') {
          U.prime.mean <- U + lmd^2/2*GU
        }
        if (Metropolis == 'Newton') {
          stop('Newton MH is only for the ordinal case!')
        }
        U.prime <- U.prime.mean + lmd * matrix(rnorm(n*D), nrow = n)
        U.prime[ind.matrix.lb] <- rtmvnorm(1, U.prime.mean[ind.matrix.lb], lmd^2 * diag(n.lb), Tau.matrix[ind.matrix.lb], rep(Inf, n.lb))
        U.prime[ind.matrix.ub] <- rtmvnorm(1, U.prime.mean[ind.matrix.ub], lmd^2 * diag(n.ub), rep(-Inf, n.ub), Tau.matrix[ind.matrix.ub])   
      }
    }
    
    W.prime <- cbind(X, U.prime)
    KWW.prime <- KernelMatrix(W.prime, Theta)
    C.prime <- KWW.prime + sigma.error^2 * diag(n)
    Ci.prime <- solve(C.prime)
    if (consider.mean) {
      mu.prime <- MeanFunction(X, U.prime, coef1, coef2)
    }else {
      mu.prime <- mu  # rep(0, n), zero mean
    }
    if (Metropolis == 'Langevin') {
      GU.prime <- Gradient(U = U.prime, mu = mu.prime, K = KWW.prime, C = C.prime, Ci = Ci.prime)
    }
    if (Metropolis == 'Newton') {
      GU.prime <- Gradient(U = U.prime, mu = mu.prime, K = KWW.prime, C = C.prime, Ci = Ci.prime)
      HU.prime <- Hessian(U = U.prime, mu = mu.prime, K = KWW.prime, C = C.prime, Ci = Ci.prime) - diag(n)
      HUi.prime <- solve(HU.prime)
      U.Newton.prime <- U.prime - lmd^2/2*HUi.prime%*%GU.prime
    }
    
    # Accept-Reject for U.prime ----
    if (Metropolis == 'Independent') {
      logp <- min(0, LogDensityRatio(y, mu.prime, C.prime, mu, C, Ci.prime, Ci))
    } else if (Metropolis == 'Random Walk') {
      logp <- min(0, LogDensityRatio(y, mu.prime, C.prime, mu, C, Ci.prime, Ci) + sum(U^2 - U.prime^2)/2)
    } else if (Metropolis == 'Langevin') {
      logp <- min(0, 
                  LogDensityRatio(y, mu.prime, C.prime, mu, C, Ci.prime, Ci) 
                  + sum(U^2 - U.prime^2)/2 
                  + (sum((U.prime - U - lmd^2/2*GU)^2) - sum((U - U.prime - lmd^2/2*GU.prime)^2))/(2*lmd^2)
                  )
    } else if (Metropolis == 'Newton') {
      logp <- min(0, 
                  LogDensityRatio(y, mu.prime, C.prime, mu, C, Ci.prime, Ci) 
                  + sum(U^2 - U.prime^2)/2
                  - (dmvnorm(U.prime, U.Newton, -lmd^2*HUi, log = T) - dmvnorm(U, U.Newton.prime, -lmd^2*HUi.prime, log = T)) 
                  )
    }
    
    if(log(runif(1)) <= logp){
      U <- U.prime
      W <- W.prime
      KWW <- KWW.prime
      C <- C.prime
      Ci <- Ci.prime
      mu <- mu.prime
      if (Metropolis == 'Langevin') {
        GU <- GU.prime
      }
      if (Metropolis == 'Newton') {
        GU <- GU.prime
        HU <- HU.prime
        HUi <- HUi.prime
        U.Newton <- U.Newton.prime
      }
      accep[k] <- 1
    }  # Otherwise, we do not update U or W.
    
    U.hist[, ((k-1)*D+1):(k*D)] <- U
    
    # Simulating Tau ----
    if (category.type == 'ordinal') {
      for (idx in 1:sum(ind)) {
        j <- which(ind)[idx]
        j.next <- which(ind)[idx+1]
        # {c = j} and {c = j.next} are two consecutive nonempty sets
        
        lb <- max(U[c == j], -Inf)
        ub <- min(U[c == j.next], Inf)
        Tau[j:(j.next-1)] <- sort(rtmvnorm(j.next-j, 0, 1, lb, ub))
        
        if (j.next == m)  # break once Tau[m-1] has been generated
          break
      }
    }else{ # category.type == 'nominal'
      for (j in 1:(m-1)) {
        Tau[j] <- rtmvnorm(1, 0, 1, max(U[which(c==j), j], -Inf), min(U[which(c>j), j], Inf))
      }
    }
    Tau.hist[, k] <- Tau
    
  }
  
  return(list(U.hist = U.hist, Tau.hist = Tau.hist, accep = accep))
}


GibbsTau <- function(seed, Tau0, U, K,
                     n, D, m, 
                     c, X, y, 
                     Theta, sigma.error) {
  # Runs a a Gibbs sampler for the posterior distribution of Tau fixing U. 
  #
  # Args:
  #   seed: random seed
  #   Tau0: initial Tau
  #   U: fixed U
  #   K: number of iterations
  #
  # Returns:
  #   MCMC samples of Tau. 
  
  set.seed(seed)
  Tau.hist <- matrix(0, nrow = m-1, ncol = K)
  
  # initializing U, W, and Tau
  W <- cbind(X, U)
  Tau <- Tau0  
  
  # Creating some index matrices for subsetting
  if (category.type == 'ordinal') {
    ind <- rep(T, m)
    for (j in 1:m){
      ind[j] <- any(c == j)
    }
    ind[1] <- T  
    ind[m] <- T  #For convenience, we can always augment u[c==1] with -Inf and u[c==m] with Inf, such that ind[1] and ind[m] are both T.
  }
  
  # Gibbs
  for(k in 1:K) {
    #### Simulating Tau ####
    if (category.type == 'ordinal') {
      for (idx in 1:sum(ind)) {
        j <- which(ind)[idx]
        j.next <- which(ind)[idx+1]
        # {c = j} and {c = j.next} are two consecutive nonempty sets
        
        lb <- max(U[c == j], -Inf)
        ub <- min(U[c == j.next], Inf)
        Tau[j:(j.next-1)] <- sort(rtmvnorm(j.next-j, 0, 1, lb, ub))
        
        if (j.next == m)  # break once Tau[m-1] has been generated
          break
      }
    }else{ # category.type == 'nominal'
      for (j in 1:(m-1)) {
        Tau[j] <- rtmvnorm(1, 0, 1, max(U[which(c==j), j], -Inf), min(U[which(c>j), j], Inf))
      }
    }
    Tau.hist[, k] <- Tau
    
  }
  
  return(Tau.hist)
}

GibbsU <- function(seed, Tau, U0, K, Metropolis, lmd,
                     n, D, m, 
                     c, X, y, 
                     Theta, sigma.error) {
  # Runs a a Gibbs sampler for the posterior distribution of U by fixing Tau. 
  #
  # Args:
  #   seed: random seed
  #   Tau: fixed Tau
  #   U0: initial U
  #   K: number of iterations
  #   Metropolis: 'Independent' for the independent proposal or 'Random Walk' for the random walk proposal
  #   lmd: step size (std of the Gaussian proposal)
  #
  # Returns:
  #   A list of MCMC samples of U, and a binary vector recording acceptances of the Metropolis steps for U. 
  
  set.seed(seed)
  U.hist <- matrix(0, nrow = n, ncol = D*K)
  accep <- rep(0, K)
  
  # initializing U, W, and Tau
  U <- U0
  W <- cbind(X, U)
  
  # Creating some index matrices for subsetting
  if (category.type == 'ordinal') {
    ind <- rep(T, m)
    for (j in 1:m){
      ind[j] <- any(c == j)
    }
    ind[1] <- T  
    ind[m] <- T  #For convenience, we can always augment u[c==1] with -Inf and u[c==m] with Inf, such that ind[1] and ind[m] are both T.
  }else{
    ind.matrix <- matrix(rep(1:D, n), nrow = n, byrow = T)
    ind.matrix.c <- matrix(rep(c, D), ncol = D)
    ind.matrix.lb <- ind.matrix - ind.matrix.c < 0  # logical matrix indicating the indices of U with lower bound constraints
    ind.matrix.ub <- ind.matrix - ind.matrix.c == 0  # logical matrix indicating the indices of U with upper bound constraints (at most 1 True in each row)
    n.lb <- sum(ind.matrix.lb)
    n.ub <- sum(ind.matrix.ub)
  }
  
  # Gibbs
  for(k in 1:K) {
    #### Simulating U ####
    # Generating the proposal U.prime
    if (category.type == 'ordinal') {
      Tau.aug <- c(-Inf, Tau, Inf)
      if (Metropolis == 'Independent') {
        U.prime <- rtmvnorm(1, rep(0,n), diag(n), Tau.aug[c], Tau.aug[c+1])
      }else if (Metropolis == 'Random Walk') {
        U.prime <- rtmvnorm(1, U, lmd^2 * diag(n), Tau.aug[c], Tau.aug[c+1])
      }else{
        stop('Unsupported Metrpolis!')
      }
      
      
    }else{  # category.type == 'nominal'
      Tau.matrix <- matrix(rep(Tau, n), nrow = n, byrow = T)
      if (Metropolis == 'Independent') {
        U.prime <- matrix(rnorm(n*D), nrow = n)
        U.prime[ind.matrix.lb] <- rtmvnorm(1, rep(0, n.lb), diag(n.lb), Tau.matrix[ind.matrix.lb], rep(Inf, n.lb))
        U.prime[ind.matrix.ub] <- rtmvnorm(1, rep(0, n.ub), diag(n.ub), rep(-Inf, n.ub), Tau.matrix[ind.matrix.ub])
      } else if (Metropolis == 'Random Walk') {
        U.prime <- matrix(U + lmd * rnorm(n*D), nrow = n)
        U.prime[ind.matrix.lb] <- rtmvnorm(1, U[ind.matrix.lb], lmd^2 * diag(n.lb), Tau.matrix[ind.matrix.lb], rep(Inf, n.lb))
        U.prime[ind.matrix.ub] <- rtmvnorm(1, U[ind.matrix.ub], lmd^2 * diag(n.ub), rep(-Inf, n.ub), Tau.matrix[ind.matrix.ub])
      }
      
    }
    
    # Accept-Reject for U
    W.prime <- cbind(X, U.prime)
    if (Metropolis == 'Independent') {
      #logp <- min(0, log(dmvnorm(y, rep(0, n), KernelMatrix(W.prime, Theta) + sigma.error^2 * diag(n))) - log(dmvnorm(y, rep(0, n), KernelMatrix(W, Theta) + sigma.error^2 * diag(n))))
      logp <- min(0, LogDensityRatio(y, rep(0, n), KernelMatrix(W.prime, Theta) + sigma.error^2 * diag(n), rep(0, n), KernelMatrix(W, Theta) + sigma.error^2 * diag(n)))
    } else if (Metropolis == 'Random Walk'){
      #logp <- min(0, log(dmvnorm(y, rep(0, n), KernelMatrix(W.prime, Theta) + sigma.error^2 * diag(n))) - log(dmvnorm(y, rep(0, n), KernelMatrix(W, Theta) + sigma.error^2 * diag(n))) + sum(U^2 - U.prime^2)/2)
      logp <- min(0, LogDensityRatio(y, rep(0, n), KernelMatrix(W.prime, Theta) + sigma.error^2 * diag(n), rep(0, n), KernelMatrix(W, Theta) + sigma.error^2 * diag(n)) + sum(U^2 - U.prime^2)/2)
    }
    
    
    if(log(runif(1)) <= logp){
      U <- U.prime
      W <- W.prime
      accep[k] <- 1
    }  # Otherwise, we do not update U or W.
    
    U.hist[, ((k-1)*D+1):(k*D)] <- U
  }
  
  return(list(U.hist = U.hist, accep = accep))
}

