#--------------------  FUNCTIONS --------------------------
KernelMatrix <- function (W, Theta, V = NULL, eps = 0) {
  # Compute the kernel matrix K_{WV} with the inverse exponential squared kernel.
  #
  # Args:
  #   W: The data matrix with each row being an observation.
  #   Theta: The last one is the overall variance, and the rest are lengthscales.
  #   V: The second data matrix, by default is the same as W.
  #
  # Returns:
  #   A kernel matrix.
  
  if (length(Theta) != ncol(W)+1) {
    stop("Wrong number of kernel parameters!")
  }
  if (any(Theta < 0)) {
    stop('Negative kernel parameter!')
  }
  
  if (is.null(V)) 
    V <- W
  if (nrow(W) <= nrow(V)) {
    K <- t(sapply(1:nrow(W), function(idx) exp(-apply((t(V)-W[idx,])^2*Theta[1:ncol(W)], 2, sum))))
  } else {
    K <- sapply(1:nrow(V), function(idx) exp(-apply((t(W)-V[idx,])^2*Theta[1:ncol(W)], 2, sum)))
  }
  K <- Theta[ncol(W)+1]*K  # Scale the kernel matrix
  if (is.null(V) && eps>0) {
    diag(K) <- diag(K) + eps # add a small perturbation to the kernel matrix for numerical inversion
  }
  return(K)
}

EstimateEvidence <- function (theta, rho,
                              U.hist,
                              IG.a.error, IG.b.error,
                              X, y, d = 1) {
  # Apply to both ordinal and nominal cases.
  # Given hyper-parameter (rho, theta), estimate the marginal likelihood up to a constant factor by Monte Carlo.
  #
  # Arg:
  #   theta: lengthscales of x and u
  #   rho: signal-to-noise ratio. That is, the GP variance is rho^2*sigma.error2.
  #   U.hist: Samples of U from p(U|X, c); d: the dimension of u
  #   IG.a.error and IG.b.error: parameters of the IG prior for sigma.error2
  #
  # Return:
  #  A vector storing the samples of marginal likelihood.
  
  n <- nrow(U.hist)
  n.mc <- ncol(U.hist)%/%d
  mar.likelihood.hist <- rep(0, n.mc)
  for (k in 1:n.mc) {
    U <- U.hist[, ((k-1)*d+1):(k*d), drop = FALSE]
    C <- KernelMatrix(W = cbind(X, U), Theta = c(theta, rho^2)) + diag(n)
    C.chol <- t(chol(C)) # the (lower-triangular) Cholesky factor of C, i.e., C = C.chol%*%t(C.chol)
    det.C <- prod(diag(C.chol))^2
    #det.C <- c(determinant(C, logarithm = FALSE)$modulus) # C>0
    r <- forwardsolve(C.chol, y)
    a.tilde <- IG.a.error + n/2
    b.tilde <- IG.b.error + sum(r^2)/2
    #b.tilde <- IG.b.error + t(y)%*%solve(C)%*%y/2
    mar.likelihood.hist[k] <- b.tilde^(-a.tilde)/sqrt(det.C)
  }
  return(mar.likelihood.hist)
}  

Predict.y <- function (d, y, X, X.star, 
                       sigma.error2.hist, U.hist, U.star.hist,
                       theta, rho, single.ustar = FALSE, separate = FALSE, seed = 888) {
  # Predict y by kriging (given U.star).
  # single.ustar: whether the goal is to predict at a fixed u.star.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Constants ----
  n.mc <- length(sigma.error2.hist)  # number of repeated samples to draw
  n.star <- nrow(X.star)  # number of new samples
  
  # Empty matrices ----
  f.predict.hist <- matrix(0, nrow = n.star, ncol = n.mc)
  y.predict.hist <- matrix(0, nrow = n.star, ncol = n.mc)
  
  # Sampler ----
  for (idx in 1:n.mc) {
    # Extract ----
    sigma.error2 <- sigma.error2.hist[idx]  
    U <- U.hist[, ((idx-1)*d+1):(idx*d), drop = FALSE]
    if (single.ustar) {
      U.star <- U.star.hist
    } else {
      U.star <- U.star.hist[, ((idx-1)*d+1):(idx*d), drop = FALSE]
    }
    
    # Draw y.star ----
    Theta <- c(theta, rho^2*sigma.error2)
    W <- cbind(X, U)
    W.star <- cbind(X.star, U.star)
    KWW <- KernelMatrix(W = W, Theta = Theta)
    C <- KWW + sigma.error2*diag(n)
    Ci <- solve(C)
    KstarW <- KernelMatrix(W = W.star, Theta = Theta, V = W)
    
    if (separate) {
      sdv <-  sqrt(rho^2*sigma.error2 - apply(KstarW, 1, function(x) t(x)%*%Ci%*%x))
      f.predict.hist[, idx] <- KstarW%*%Ci%*%y + rnorm(n.star)*sdv
    } else {
      Kstar <- KernelMatrix(W = W.star, Theta = Theta)
      var.star <- Kstar - KstarW%*%Ci%*%t(KstarW)
      var.star <- as.symmetric.matrix(var.star)  #in case of numerical error
      f.predict.hist[, idx] <- rmvnorm(1, KstarW%*%Ci%*%y, var.star)  # kriging
    }
    y.predict.hist[, idx] <- f.predict.hist[, idx] + rnorm(n.star, 0, sqrt(sigma.error2))
  }
  
  # Return ----
  return(list(f.predict.hist = f.predict.hist,
              y.predict.hist = y.predict.hist))
}

EnergyDistance <- function (X, Y) {
  # X: predictions (a vector); Y: samples from the target
  # 2*E|X-Y| - E|X-X'| - E|Y-Y'|
  n <- length(X)
  m <- length(Y)
  sxy <- 0
  sxx <- 0
  syy <- 0
  for (i in 1:n) {
    for (j in 1:m) {
      sxy <- sxy + abs(X[i] - Y[j])
    }
  }
  for (i in 1:n) {
    for (j in 1:n) {
      sxx <- sxx + abs(X[i] - X[j])
    }
  }
  for (i in 1:m) {
    for (j in 1:m) {
      syy <- syy + abs(Y[i] - Y[j])
    }
  }
  return(2*sxy/(n*m) - sxx/n^2 - syy/m^2)
}

#---------------------- Others -----------------------------
GetMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}