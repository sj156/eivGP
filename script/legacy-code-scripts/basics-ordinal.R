#------------- Packages required ------------
# library('GA') for visualization
# library("mvtnorm") for simulating multivariate normal

#--------------------  SIMULATION --------------------------
ResponseFunction <- function (X, U, code) {
  if (code == 'sin') {
    return(sin(U))
  } else if (code == 'cos') {
    return(cos(U))
  } else if (code == 'quad') {
    return((X + U)^2)
  } else if (code == 'doppler') {
    return(sqrt(U*(1-U))*sin(2.1*pi/(U+0.05)))
  } else if (code == 'sin-step') {
    U.ref <- c(0:4)*pi/2
    U <- apply(U, 1, function(u) U.ref[which.min(abs(u - U.ref))])
    return(sin(U))
  } else if (code == 'quad-step') {
    U.ref <- c(-1.5, -1, 0, 1, 1.5)
    U <- apply(U, 1, function(u) U.ref[which.min(abs(u - U.ref))])
    return((X + U)^2)
  } else {
    stop('Invalid code!')
  }
}

PlotResponse <- function (code, nx = 100, lb.u = -2, ub.u = 2, lb.x = -2, ub.x = 2) {
  # Draw a perspective plot for a given function f(x,u). Only univariate x and u are allowed. 
  x <- seq(lb.x, ub.x, length = nx)
  u <- seq(lb.u, ub.u, length = nx)
  W <- expand.grid(x, u)
  y <- ResponseFunction(matrix(W[, 1], ncol = 1), matrix(W[, 2], ncol = 1), code)
  persp3D(x, u, matrix(y, ncol = nx), xlab = "x", ylab = "u", zlab = "y")
}

#--------------------  INFERENCE --------------------------
# For vectors, at least one of the arg should be a vector
DrawTrunNormal <- function (mu, sigma, lb, ub) {
  # Simulate independent truncated normals by inverse probability integral transform.
  # 
  # Args:
  #   mu: vector of means; sigma: vector of std; lb: vector of lower bounds; ub: vector of upper bounds
  
  n <- max(length(mu), length(sigma), length(lb), length(ub))
  u <- runif(n, pnorm(lb, mu, sigma), pnorm(ub, mu, sigma))
  return(qnorm(u, mu, sigma))
}

Gradient <- function (U, r, phi, K, length.scale.u, X.tilde, beta, sigma.u2) {
  # Calculate the gradient of the log conditional density w.r.t U 
  
  C.tilde <- (phi%*%t(phi) - r%*%t(r))*K
  D.U <- 2*length.scale.u*(diag(apply(C.tilde, 1, sum)) - C.tilde)%*%U - (U - X.tilde%*%beta)/sigma.u2
  return(D.U)
}

Gibbs <- function (seed, num.iter, lmd, 
                   beta0, sigma.u20, sigma.error20, Tau0, U0, 
                   X, c, y, u.obs.idx, U.obs, u.lb = -Inf, u.ub = Inf,
                   v0, g, IG.a.u, IG.b.u, IG.a.error, IG.b.error,
                   theta, rho, simu.f = TRUE, design.map = FALSE) {
  
  if(!is.null(seed))
    set.seed(seed)
  
  # Compute constants ----
  n <- length(c)
  m <- length(Tau0) + 1
  d.x <- ncol(X)
  n.obs <- length(u.obs.idx)
  if (n.obs > 0) {
    u.mis.idx <- c(1:n)[-u.obs.idx] # the index for missing u
  } else {
    u.mis.idx <- c(1:n)
  }
  ##
  X.tilde <- cbind(1, X)
  XX.tilde <- t(X.tilde)%*%X.tilde
  V0 <- matrix(0, d.x+1, d.x+1)
  V0[1,1] <- v0
  V0[-1,-1] <- g*XX.tilde[-1,-1]
  Vn <- V0 + XX.tilde 
  Vn.inv <- solve(Vn)
  ## for the MAP estimation of u~X
  if (design.map) {
    X.tilde.obs <- X.tilde[u.obs.idx, ]
    XX.tilde.obs <- t(X.tilde.obs)%*%X.tilde.obs
    V0.obs <- matrix(0, d.x+1, d.x+1)
    V0.obs[1,1] <- v0
    V0.obs[-1,-1] <- g*XX.tilde.obs[-1,-1]
    Vn.obs <- V0.obs + XX.tilde.obs
    Vn.obs.inv <- solve(Vn.obs)
    mu.n.obs <- Vn.obs.inv%*%t(X.tilde.obs)%*%U.obs
    beta.map <- mu.n.obs
    sigma.u2.map <- (IG.b.u + (t(U.obs)%*%U.obs - t(mu.n.obs)%*%Vn.obs%*%mu.n.obs)/2)/(IG.a.u + n.obs/2 + (d.x+1)/2 + 1)
  }
  
  # Create matrices to store Gibbs samples ----
  beta.hist <- matrix(0, nrow = d.x + 1, ncol = num.iter)
  sigma.u2.hist <- rep(0, num.iter)
  sigma.error2.hist <- rep(0, num.iter)
  Tau.hist <- matrix(0, nrow = m-1, ncol = num.iter)
  U.hist <- matrix(0, nrow = n, ncol = num.iter)
  accep <- rep(0, num.iter)
  f.hist <- matrix(0, nrow = n, ncol = num.iter)
  
  # Initialize ----
  beta <- beta0
  sigma.u2 <- sigma.u20
  sigma.error2 <- sigma.error20
  Tau <- Tau0  
  U <- U0
  if (n.obs > 0) {
    U[u.obs.idx] <- U.obs  # set observed u
  }
  
  # Variables kept latest ----
  Theta <- c(theta, rho^2*sigma.error2)
  W <- cbind(X, U)
  KWW <- KernelMatrix(W, Theta)
  C <- KWW + sigma.error2*diag(n)
  C.chol <- t(chol(C)) # the (lower-triangular) Cholesky factor of C, i.e., C = C.chol%*%t(C.chol)
  r <- backsolve(t(C.chol), forwardsolve(C.chol, y)) # i.e., r = C^{-1}y
  
  # Gibbs ----
  for (k in 1:num.iter) {
    # Simulate {beta, sigma.u2} ----
    if (design.map) {
      sigma.u2 <- sigma.u2.map
      beta <- beta.map
    } else {
      mu.n <- Vn.inv%*%t(X.tilde)%*%U
      sigma.u2 <- 1/rgamma(1, shape = IG.a.u + n/2, rate = IG.b.u + (t(U)%*%U - t(mu.n)%*%Vn%*%mu.n)/2)
      beta <- c(rmvnorm(1, mu.n, sigma.u2*Vn.inv))
    }
    ## 
    sigma.u2.hist[k] <- sigma.u2
    beta.hist[, k] <- beta 
    
    # Simulate sigma.error2 ----
    sigma.error2.prime <- 1/rgamma(1, shape = IG.a.error + n/2, 
                                   rate = IG.b.error + sigma.error2/2*t(y)%*%r) # buffering
    
    # Maintain variables
    Theta <- c(theta, rho^2*sigma.error2.prime)
    KWW <- (sigma.error2.prime/sigma.error2)*KWW
    C <- (sigma.error2.prime/sigma.error2)*C
    C.chol <- sqrt(sigma.error2.prime/sigma.error2)*C.chol
    r <- (sigma.error2/sigma.error2.prime)*r
    sigma.error2 <- sigma.error2.prime  # remove buffering
    sigma.error2.hist[k] <- sigma.error2
    
    # Simulate Tau ----
    for (j in 1:(m-1)) {
      # Tau[j] <- rmvnorm(1, 0, 1, max(U[c == j]), min(U[c == j+1]))
      Tau[j] <- runif(1, max(U[c == j]), min(U[c == j+1]))
    }
    Tau.hist[, k] <- Tau
    
    # Simulate phi and LMC ----
    phi <- backsolve(t(C.chol), rnorm(n))
    GU <- Gradient(U = U, r = r, phi = phi, K = KWW, length.scale.u = theta[d.x + 1], X.tilde = X.tilde, beta = beta, sigma.u2 = sigma.u2)
    Tau.aug <- c(u.lb, Tau, u.ub)
    U.mis <- U[u.mis.idx]
    U.mis.prime <- DrawTrunNormal(U.mis + lmd^2/2*GU[u.mis.idx], rep(lmd, length(u.mis.idx)), Tau.aug[c[u.mis.idx]], Tau.aug[c[u.mis.idx]+1])
    
    # Maintain variables ----
    U.prime <- U  # note that U[u.obs.idx] has already been set to U.obs when initialized
    U.prime[u.mis.idx] <- U.mis.prime  # set U.prime
    W.prime <- cbind(X, U.prime)
    KWW.prime <- KernelMatrix(W.prime, Theta)
    C.prime <- KWW.prime + sigma.error2*diag(n)
    C.prime.chol <- t(chol(C.prime))
    r.prime <- backsolve(t(C.prime.chol), forwardsolve(C.prime.chol, y))
    GU.prime <- Gradient(U = U.prime, r = r.prime, phi = phi, K = KWW.prime, length.scale.u = theta[d.x + 1],
                         X.tilde = X.tilde, beta = beta, sigma.u2 = sigma.u2) 
    
    # Calculate the log acceptance probability ----
    p.u <- pnorm(Tau.aug[c[u.mis.idx]+1], mean = U.mis + lmd^2/2*GU[u.mis.idx], sd = lmd) - pnorm(Tau.aug[c[u.mis.idx]], mean = U.mis + lmd^2/2*GU[u.mis.idx], sd = lmd) 
    p.u.prime <- pnorm(Tau.aug[c[u.mis.idx]+1], mean = U.mis.prime + lmd^2/2*GU.prime[u.mis.idx], sd = lmd) - pnorm(Tau.aug[c[u.mis.idx]], mean = U.mis.prime + lmd^2/2*GU.prime[u.mis.idx], sd = lmd) 
    log.ratio <- try(sum(log(p.u)) - sum(log(p.u.prime)), TRUE)
    if (inherits(log.ratio, "try-error")) {
      log.ratio <- 0
    }
    logp <- min(0, -(t(y)%*%(r.prime - r) + t(phi)%*%(C.prime - C)%*%phi)/2 - 
                  sum((U.prime - X.tilde%*%beta)^2 - (U - X.tilde%*%beta)^2)/(2*sigma.u2) + 
                  sum((U.mis.prime - U.mis - lmd^2/2*GU[u.mis.idx])^2 - (U.mis - U.mis.prime - lmd^2/2*GU.prime[u.mis.idx])^2)/(2*lmd^2) +
                  log.ratio)
    
    # Accept-Reject for U.prime ----
    if(log(runif(1)) <= logp){  # Accept U.prime.
      U <- U.prime
      W <- W.prime
      KWW <- KWW.prime
      C <- C.prime
      C.chol <- C.prime.chol
      r <- r.prime
      accep[k] <- 1
    }  # Otherwise, do not update U.
    
    U.hist[, k] <- U
    # Simulate f (costly) ----
    if (simu.f) {
      LiK <- apply(KWW, 2, function(b) forwardsolve(C.chol, b))  # i.e., solve(C.chol)%*%KWW
      f.hist[, k] <- rmvnorm(1, KWW%*%r, KWW - t(LiK)%*%LiK)
    } else {
      f.hist[, k] <- rep(0, n)
    }
  }
  
  # return ----
  return(list(beta.hist = beta.hist,
              sigma.u2.hist = sigma.u2.hist,
              sigma.error2.hist = sigma.error2.hist,
              Tau.hist = Tau.hist,
              U.hist = U.hist,
              accep = accep,
              f.hist = f.hist))
}

Gibbs.EB <- function (seed, n.mc, 
                      beta0, sigma.u20, Tau0, U0, 
                      X, c, u.obs.idx, U.obs, u.lb = -Inf, u.ub = Inf,
                      v0, g, IG.a.u, IG.b.u) {
  # Gibbs sampler for {U.mis, Tau, beta, sigma.u2} given {U.obs, X, c}.
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Compute constants----
  n <- nrow(X)
  m <- length(Tau0) + 1
  d.x <- ncol(X)
  X.tilde <- cbind(1, X)
  XX.tilde <- t(X.tilde)%*%X.tilde
  V0 <- matrix(0, d.x+1, d.x+1)
  V0[1,1] <- v0
  V0[-1,-1] <- g*XX.tilde[-1,-1]
  Vn <- V0 + XX.tilde 
  Vn.inv <- solve(Vn)
  if (length(u.obs.idx) > 0) {
    u.mis.idx <- c(1:n)[-u.obs.idx] # the index for missing u
  } else {
    u.mis.idx <- c(1:n)
  }
  
  # Initialize ----
  beta <- beta0
  sigma.u2 <- sigma.u20
  Tau <- Tau0
  U <- U0
  if (length(u.obs.idx) > 0) {
    U[u.obs.idx] <- U.obs  # set observed u
  }
  
  beta.hist <- matrix(0, nrow = d.x+1, ncol = n.mc)
  sigma.u2.hist <- rep(0, n.mc)
  Tau.hist <- matrix(0, nrow = m-1, ncol = n.mc)
  U.hist <- matrix(0, nrow = n, ncol = n.mc)
  
  # Gibbs ----
  for (k in 1:n.mc) {
    # Simulate beta & sigma.u2 (SAME as inference) ----
    mu.n <- Vn.inv%*%t(X.tilde)%*%U
    sigma.u2 <- 1/rgamma(1, shape = IG.a.u + n/2, rate = IG.b.u + (t(U)%*%U - t(mu.n)%*%Vn%*%mu.n)/2)
    beta <- c(rmvnorm(1, mu.n, sigma.u2*Vn.inv))
    sigma.u2.hist[k] <- sigma.u2
    beta.hist[, k] <- beta 
    
    # Simulate Tau (SAME as inference) ----
    for (j in 1:(m-1)) {
      # Tau[j] <- rmvnorm(1, 0, 1, max(U[c == j]), min(U[c == j+1]))
      Tau[j] <- runif(1, max(U[c == j]), min(U[c == j+1]))
    }
    Tau.hist[, k] <- Tau
    
    # Sample (missing) U ----
    Tau.aug <- c(u.lb, Tau, u.ub)
    U[u.mis.idx] <- DrawTrunNormal(X.tilde[u.mis.idx, ]%*%beta, rep(sqrt(sigma.u2), length(u.mis.idx)), Tau.aug[c[u.mis.idx]], Tau.aug[c[u.mis.idx]+1])  # U[u.obs.idx] has already been set to U.obs when initialized
    U.hist[, k] <- U
  }
  
  # Return ----
  return(list(beta.hist = beta.hist, 
              sigma.u2.hist = sigma.u2.hist,
              Tau.hist = Tau.hist,
              U.hist = U.hist))
}

Gibbs.U.star <- function (X.star, c.star, beta.hist, sigma.u2.hist, Tau.hist, 
                          u.lb = -Inf, u.ub = Inf, seed = 222) {
  # Simulate u* given {X*, c*, beta, sigma.u2, Tau}
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Constants ----
  n.mc <- ncol(Tau.hist)  # number of repeated samples to draw
  n.star <- length(c.star)  # number of new samples
  X.star.tilde <- cbind(1, X.star)
  
  # Empty matrices ----
  U.star.hist <- matrix(0, nrow = n.star, ncol = n.mc)

  # Sampler ----
  for (idx in 1:n.mc) {
    # extract beta, sigma.u2, Tau ----
    beta <- beta.hist[, idx]
    sigma.u2 <- sigma.u2.hist[idx]
    Tau <- Tau.hist[, idx]
    Tau.aug <- c(u.lb, Tau, u.ub)
    
    # draw U.star ----
    U.star <- DrawTrunNormal(X.star.tilde%*%beta, rep(sqrt(sigma.u2), n.star), Tau.aug[c.star], Tau.aug[c.star+1])
    U.star.hist[, idx] <- U.star
  }
  
  # Return ----
  return(U.star.hist = U.star.hist)
}

Predict.y.oracle <- function (code, X.star, c.star, 
                              beta, sigma.u2, sigma.error2, Tau, u.lb, u.ub,
                              uniform = FALSE, n.mc = 10000, seed = 999) {
  # Predict (generate) the output y (and f) using the true parameters (assuming the DGP is true)
  # uniform: if true, the design of u is a uniform distribution independent of x. Otherwise, it follows DGP (i.e., Gaussian).
  
  # Initialize ----
  set.seed(seed)
  n.star <- length(c.star)
  U.star.hist <- matrix(0, nrow = n.star, ncol = n.mc)
  f.predict.hist <- matrix(0, nrow = n.star, ncol = n.mc)
  y.predict.hist <- matrix(0, nrow = n.star, ncol = n.mc)
  
  # Constants ----
  X.star.tilde <- cbind(1, X.star)
  
  # Predict ----
  for (idx in 1:n.mc) {
    # Draw U.star ----
    Tau.aug <- c(u.lb, Tau, u.ub)
    if (uniform) {
      U.star <- runif(n.star, Tau.aug[c.star], Tau.aug[c.star+1])
    } else {
      U.star <- DrawTrunNormal(X.star.tilde%*%beta, rep(sqrt(sigma.u2), n.star), Tau.aug[c.star], Tau.aug[c.star+1])
    }
    U.star.hist[, idx] <- U.star
    
    # Draw y.star and f.star ----
    f.predict.hist[, idx] <- ResponseFunction(X.star, as.matrix(U.star), code) 
    y.predict.hist[, idx] <- f.predict.hist[, idx] + rnorm(n.star, 0, sqrt(sigma.error2))
  }
  
  # Return ----
  return(list(U.star.hist = U.star.hist,
              y.predict.hist = y.predict.hist, 
              f.predict.hist = f.predict.hist))
}

GetInitial <- function (n, m, c, u.obs.idx, U.obs, u.lb = -Inf, u.ub = Inf, seed = 555) {
  # Generate initial value of U and Tau for Gibbs.
  # When U.obs is given, it is allowed that u are totally missing for some classes.
  U <- rep(0, n)
  if (length(u.obs.idx) > 0) {
    Tau.aug <- c(u.lb, rep(0, m-1), u.ub)
    for (j in 1:(m-1)) {
      Tau.aug[j+1] <- DrawTrunNormal(0, 1, max(Tau.aug[j], U.obs[c[u.obs.idx] <= j]), min(u.ub, U.obs[c[u.obs.idx] > j]))  # tau_j
    }
    U <- DrawTrunNormal(0, 1, Tau.aug[c], Tau.aug[c+1])
    U[u.obs.idx] <- U.obs
    Tau <- Tau.aug[2:m]
  } else {
    Tau <- sort(rnorm(m-1, u.lb, u.ub))
    Tau.aug <- c(u.lb, Tau, u.ub)
    U <- DrawTrunNormal(0, 1, Tau.aug[c], Tau.aug[c+1])
  }
  return(list(U0 = U, Tau0 = Tau, Tau.aug0 = Tau.aug))
}

