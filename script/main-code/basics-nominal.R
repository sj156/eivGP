#------------- PACKAGES REQUIRED ------------
# library("pgdraw") for simulating polya-Gamma
# library("mvtnorm") for simulating multivariate normal
# library('LaplacesDemon') for sampling inverse Wishart and  Matrix Normal

#------------- INFERENCE ---------------
Gibbs <- function (seed, num.iter, lmd,
                   B.tilde0, Sigma.u0, sigma.error20, gamma.mat0, U0, 
                   m, c, X, y, u.obs.idx, U.obs,
                   v0, g, nu0, Psi0, q0, q.gamma, IG.a.error, IG.b.error,
                   theta, rho, simu.f = TRUE, design.map = FALSE) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Compute constants ----
  n <- nrow(X)
  d.x <- ncol(X)
  d <- ncol(U0)
  n.obs <- length(u.obs.idx)
  if (n.obs > 0) {
    u.mis.idx <- c(1:n)[-u.obs.idx] # the index for missing u
  } else {
    u.mis.idx <- c(1:n)
  }
  ##
  bin.encoding <- BinaryEncode(c, m)
  N.mat <- bin.encoding$N.mat
  kappa.mat <- bin.encoding$kappa.mat
  ##
  X.tilde <- cbind(1, X)
  XX.tilde <- t(X.tilde)%*%X.tilde
  V0 <- matrix(0, d.x+1, d.x+1)
  V0[1,1] <- v0
  V0[-1,-1] <- g*XX.tilde[-1,-1]
  Vn <- V0 + XX.tilde
  Vn.inv <- as.symmetric.matrix(solve(Vn)) 
  Q0.gamma <- diag(c(q0, rep(q.gamma, d)))
  ## for the MAP estimation of U~X
  if (design.map) {
    nun.obs <- nu0 + n.obs
    X.tilde.obs <- X.tilde[u.obs.idx, ]
    XX.tilde.obs <- t(X.tilde.obs)%*%X.tilde.obs
    V0.obs <- matrix(0, d.x+1, d.x+1)
    V0.obs[1,1] <- v0
    V0.obs[-1,-1] <- g*XX.tilde.obs[-1,-1]
    Vn.obs <- V0.obs + XX.tilde.obs
    Vn.obs.inv <- solve(Vn.obs)
    Mn.obs <- Vn.obs.inv%*%t(X.tilde.obs)%*%U.obs
    Psin.obs <- Psi0 + t(U.obs)%*%U.obs - t(Mn.obs)%*%Vn.obs%*%Mn.obs
    B.tilde.map <- Mn.obs
    Sigma.u.map <- Psin.obs/(nun.obs+d.x+1+d+1)
  }
  
  # Create matrices to store Gibbs samples ----
  B.tilde.hist <- matrix(0, d.x+1, d*num.iter)
  Sigma.u.hist <- matrix(0, d, d*num.iter)
  sigma.error2.hist <- rep(0, num.iter)
  gamma.mat.hist <- matrix(0, m-1, (d+1)*num.iter)
  U.hist <- matrix(0, n, d*num.iter)
  accep <- rep(0, num.iter)
  f.hist <- matrix(0, n, num.iter)
  
  # Initialize ----
  B.tilde <- B.tilde0
  Sigma.u <- Sigma.u0
  sigma.error2 <- sigma.error20
  gamma.mat <- gamma.mat0
  U <- U0
  if (n.obs > 0) {
    U[u.obs.idx, ] <- U.obs  # set observed u. Only U[u.mis.idx,] will be updated.
  }
  ##
  omega.mat <- matrix(0, n, m-1)  # by default 0
  mn.u <- matrix(0, n, d)
  Qn.u <- array(0, c(d,d,n))
  Qn.u.inv <- array(0, c(d,d,n))
  
  # Variables kept latest ----
  Theta <- c(theta, rho^2*sigma.error2)
  W <- cbind(X, U)
  KWW <- KernelMatrix(W, Theta)
  C <- KWW + sigma.error2*diag(n)
  C.chol <- t(chol(C))  # the (lower-triangular) Cholesky factor of C, i.e., C = C.chol%*%t(C.chol)
  r <- backsolve(t(C.chol), forwardsolve(C.chol, y))  # i.e., r = C^{-1}%*%y
  
  # Gibbs ----
  for (k.iter in 1:num.iter) {
    # Simulate {B.tilde, Sigma.u} ----
    if (design.map) {
      Sigma.u <- Sigma.u.map
      B.tilde <- B.tilde.map
    } else {
      nun <- nu0 + n
      Mn <- Vn.inv%*%t(X.tilde)%*%U
      Psin <- Psi0 + t(U)%*%U - t(Mn)%*%Vn%*%Mn
      Sigma.u <- rinvwishart(nun, Psin)
      B.tilde <- rmatrixnorm(Mn, Vn.inv, Sigma.u)
    }
    ##
    Sigma.u.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- Sigma.u
    B.tilde.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- B.tilde
    
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
    sigma.error2.hist[k.iter] <- sigma.error2
    
    # Simulate omega ----
    U.tilde <- cbind(1, U)
    psi.mat <- U.tilde%*%t(gamma.mat)
    omega.mat[N.mat == 1] <- pgdraw(1, psi.mat[N.mat == 1])  # only update those with N_{ij}=1, others remain as the initial 0.
    
    # Simulate gamma ----
    # U.tilde <- cbind(1, U) # remove this if U.tilde has been updated
    for (j in 1:(m-1)) {
      Qn.gamma <- Q0.gamma + t(U.tilde)%*%diag(omega.mat[, j])%*%U.tilde
      Qn.gamma.inv <- solve(Qn.gamma)
      mn.gamma <- Qn.gamma.inv%*%t(U.tilde)%*%kappa.mat[, j]
      gamma.mat[j, ] <- rmvnorm(1, mn.gamma, Qn.gamma.inv)
    }
    gamma.mat.hist[, ((k.iter-1)*(d+1)+1):(k.iter*(d+1))] <- gamma.mat
    
    # Simulate phi and LMC for U ----
    phi <- backsolve(t(C.chol), rnorm(n))
    Sigma.u.inv <- solve(Sigma.u)
    for (i in 1:n) {
      Qn.u[,,i] <-  Sigma.u.inv + t(gamma.mat[,-1])%*%diag(omega.mat[i,])%*%gamma.mat[,-1]
      Qn.u.inv[,,i] <- solve(Qn.u[,,i])
      mn.u[i, ] <- Qn.u.inv[,,i]%*%(Sigma.u.inv%*%t(B.tilde)%*%X.tilde[i,] + t(gamma.mat[,-1])%*%(kappa.mat[i,] - omega.mat[i,]*gamma.mat[, 1]))
    }
    GU <- Gradient(U, phi, r, KWW, mn.u, Qn.u, theta[d.x+1])
    U.mis <- U[u.mis.idx, ]
    U.mis.prime <- U.mis + lmd^2/2*GU[u.mis.idx, ] + lmd*rnorm(length(u.mis.idx)*d)
    
    # Maintain variables ----
    U.prime <- U  # note that U[u.obs.idx, ] has already been set to U.obs when initialized
    U.prime[u.mis.idx, ] <- U.mis.prime  # update the missing part of U
    W.prime <- cbind(X, U.prime)
    KWW.prime <- KernelMatrix(W.prime, Theta)
    C.prime <- KWW.prime + sigma.error2*diag(n)
    C.prime.chol <- t(chol(C.prime))
    r.prime <- backsolve(t(C.prime.chol), forwardsolve(C.prime.chol, y))
    GU.prime <- Gradient(U.prime, phi, r.prime, KWW.prime, mn.u, Qn.u, theta[d.x+1])
    
    # Calculate the log acceptance probability ----
    logp.gau <- sapply(u.mis.idx, function(i) {
      t(U.prime[i,] - mn.u[i,])%*%Qn.u[,,i]%*%(U.prime[i,] - mn.u[i,]) - t(U[i,] - mn.u[i,])%*%Qn.u[,,i]%*%(U[i,] - mn.u[i,])
    })
    logp <- min(0, -sum(logp.gau)/2 - 
                  (t(y)%*%(r.prime - r) + t(phi)%*%(C.prime - C)%*%phi)/2 + 
                  sum((U.mis.prime - U.mis - lmd^2/2*GU[u.mis.idx,])^2 - (U.mis - U.mis.prime - lmd^2/2*GU.prime[u.mis.idx,])^2)/(2*lmd^2))
    
    # Accept-reject for U.prime ----
    if (log(runif(1)) <= logp) {  # Accept U.prime.
      U <- U.prime
      W <- W.prime
      KWW <- KWW.prime
      C <- C.prime
      C.chol <- C.prime.chol
      r <- r.prime
      accep[k.iter] <- 1
    }  # Otherwise, do not update U.
    U.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- U
    
    # Simulate f (costly) ----
    if (simu.f) {
      LiK <- apply(KWW, 2, function(b) forwardsolve(C.chol, b))  # i.e., solve(C.chol)%*%KWW
      f.hist[, k.iter] <- rmvnorm(1, KWW%*%r, KWW - t(LiK)%*%LiK)
    } else {
      f.hist[, k.iter] <- rep(0, n)
    }
  }
  
  # Return ----
  return(list(B.tilde.hist = B.tilde.hist,
              Sigma.u.hist = Sigma.u.hist,
              sigma.error2.hist = sigma.error2.hist,
              gamma.mat.hist = gamma.mat.hist,
              U.hist = U.hist, 
              accep = accep,
              f.hist = f.hist))
} 

Gradient <- function (U, phi, r, K, mn.u, Qn.u, length.scale.u) {
  n <- nrow(U)
  C.tilde <- (phi%*%t(phi) - r%*%t(r))*K
  D.U <- 2*length.scale.u*(diag(apply(C.tilde, 1, sum)) - C.tilde)%*%U  # the non-Gaussian part
  for (i in 1:n) {
    D.U[i, ] <- D.U[i, ] - Qn.u[,,i]%*%(U[i,] - mn.u[i,])  # the Gaussian part
  }
  return(D.U)
}

Gibbs.EB <- function (seed, n.mc,
                      B.tilde0, Sigma.u0, gamma.mat0, U0,
                      X, c, m, u.obs.idx, U.obs,
                      v0, g, nu0, Psi0, q0, q.gamma) {
  # Gibbs sampler for {U.mis, B.tilde, Sigma.u, gamma.mat, omega.mat} given {U.obs, X, c}
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Compute constants ----
  n <- nrow(X)
  d.x <- ncol(X)
  d <- ncol(U0)
  if (length(u.obs.idx) > 0) {
    u.mis.idx <- c(1:n)[-u.obs.idx] # the index for missing u
  } else {
    u.mis.idx <- c(1:n)
  }
  ##
  bin.encoding <- BinaryEncode(c, m)
  N.mat <- bin.encoding$N.mat
  kappa.mat <- bin.encoding$kappa.mat
  ##
  X.tilde <- cbind(1, X)
  XX.tilde <- t(X.tilde)%*%X.tilde
  V0 <- matrix(0, d.x+1, d.x+1)
  V0[1,1] <- v0
  V0[-1,-1] <- g*XX.tilde[-1,-1]
  Vn <- V0 + XX.tilde
  Vn.inv <- solve(Vn)
  Q0.gamma <- diag(c(q0, rep(q.gamma, d)))
  
  # Create matrices to store Gibbs samples ----
  B.tilde.hist <- matrix(0, d.x+1, d*n.mc)
  Sigma.u.hist <- matrix(0, d, d*n.mc)
  gamma.mat.hist <- matrix(0, m-1, (d+1)*n.mc)
  U.hist <- matrix(0, n, d*n.mc)
  
  # Initialize ----
  B.tilde <- B.tilde0
  Sigma.u <- Sigma.u0
  gamma.mat <- gamma.mat0
  U <- U0
  if (length(u.obs.idx) > 0) {
    U[u.obs.idx, ] <- U.obs  # set observed u. Only U[u.mis.idx, ] will be sampled.
  }
  ##
  omega.mat <- matrix(0, n, m-1)  # by default 0
  
  # Gibbs ----
  for (k.iter in 1:n.mc) {
    # Simulate {B.tilde, Sigma.u} (same as inference) ----
    nun <- nu0 + n
    Mn <- Vn.inv%*%t(X.tilde)%*%U
    Psin <- Psi0 + t(U)%*%U - t(Mn)%*%Vn%*%Mn
    ##
    Sigma.u <- rinvwishart(nun, Psin)
    B.tilde <- rmatrixnorm(Mn, Vn.inv, Sigma.u)
    ##
    Sigma.u.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- Sigma.u
    B.tilde.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- B.tilde 
    
    # Simulate omega (same as inference) ----
    U.tilde <- cbind(1, U)
    psi.mat <- U.tilde%*%t(gamma.mat)
    omega.mat[N.mat == 1] <- pgdraw(1, psi.mat[N.mat == 1])  # only update those with N_{ij}=1, others remain as the initial 0.
    
    # Simulate gamma (same as inference) ----
    U.tilde <- cbind(1, U) 
    for (j in 1:(m-1)) {
      Qn.gamma <- Q0.gamma + t(U.tilde)%*%diag(omega.mat[, j])%*%U.tilde
      Qn.gamma.inv <- solve(Qn.gamma)
      mn.gamma <- Qn.gamma.inv%*%t(U.tilde)%*%kappa.mat[, j]
      gamma.mat[j, ] <- rmvnorm(1, mn.gamma, Qn.gamma.inv)
    }
    gamma.mat.hist[, ((k.iter-1)*(d+1)+1):(k.iter*(d+1))] <- gamma.mat
    
    # Simulate missing U ----
    Sigma.u.inv <- solve(Sigma.u)
    for (i in u.mis.idx) {
      Qn.u <-  Sigma.u.inv + t(gamma.mat[,-1])%*%diag(omega.mat[i,])%*%gamma.mat[,-1]
      Qn.u.inv <- solve(Qn.u)
      mn.u <- Qn.u.inv%*%(Sigma.u.inv%*%t(B.tilde)%*%X.tilde[i,] + t(gamma.mat[,-1])%*%(kappa.mat[i,] - omega.mat[i,]*gamma.mat[,1]))
      U[i, ] <- rmvnorm(1, mn.u, Qn.u.inv)
    }
    U.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- U
  }
  
  # Return ----
  return(list(B.tilde.hist = B.tilde.hist,
              Sigma.u.hist = Sigma.u.hist,
              gamma.mat.hist = gamma.mat.hist,
              U.hist = U.hist))
}

Gibbs.U <- function (X, c, B.tilde.hist, Sigma.u.hist, gamma.mat.hist, 
                     U0, num.iter, seed = 222) {
  # Simulate new U* given new {X*, c*} and {B.tilde, Sigma.u, gamma.mat}. 
  # For each {B.tilde, Sigma.u, gamma.mat}, run a inner-loop Gibbs for {U*, omega*} with num.iter iterations.
  # For convenience in coding, X* and c* are just named X and c.
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Constants ----
  n <- length(c)
  m <- nrow(gamma.mat.hist) + 1
  d <- ncol(U0)
  n.mc <- ncol(Sigma.u.hist)/d   
  ##
  bin.encoding <- BinaryEncode(c, m)
  N.mat <- bin.encoding$N.mat
  kappa.mat <- bin.encoding$kappa.mat
  ##
  X.tilde <- cbind(1, X)
  
  # Create empty matrices ----
  U.hist <- matrix(0, n, d*n.mc)
  
  # Sampler ----
  for (idx in 1:n.mc) {
    # initialize ----
    U <- U0
    omega.mat <- matrix(0, n, m-1)  # initialize with 0
    
    # extract {B.tilde, Sigma.u, gamma.mat} ----
    B.tilde <- B.tilde.hist[, ((idx-1)*d+1):(idx*d)]
    Sigma.u <- Sigma.u.hist[, ((idx-1)*d+1):(idx*d)]
    Sigma.u.inv <- solve(Sigma.u)
    gamma.mat <- gamma.mat.hist[, ((idx-1)*(d+1)+1):(idx*(d+1))]
    
    # inner-loop Gibbs ----
    for (k.iter in 1:num.iter) {
      # draw omega
      U.tilde <- cbind(1, U)
      psi.mat <- U.tilde%*%t(gamma.mat)
      omega.mat[N.mat == 1] <- pgdraw(1, psi.mat[N.mat == 1])  # only update those with N_{ij}=1, others remain as the initial 0.
      
      # draw U
      for (i in 1:n) {
        Qn.u <-  Sigma.u.inv + t(gamma.mat[,-1])%*%diag(omega.mat[i,])%*%gamma.mat[,-1]
        Qn.u.inv <- solve(Qn.u)
        mn.u <- Qn.u.inv%*%(Sigma.u.inv%*%t(B.tilde)%*%X.tilde[i,] + t(gamma.mat[,-1])%*%(kappa.mat[i,] - omega.mat[i,]*gamma.mat[,1]))
        U[i, ] <- rmvnorm(1, mn.u, Qn.u.inv)
      }
    }
    
    # save the latest U ----
    U.hist[, ((idx-1)*d+1):(idx*d)] <- U
  }
  
  # Return ----
  return(U.hist)
}

Predict.U.oracle <- function (X, c, B.tilde, Sigma.u, f, uniform = FALSE, u.lb = NULL, u.ub = NULL, n.mc = 1000, seed = 999) {
  # This function predicts U using the true parameter and classifier
  # uniform: if true, the design of u is a uniform distribution independent of x. Otherwise, it follows DGP (i.e., Gaussian).
  set.seed(seed)
  
  # Initialize ----
  n <- length(c)
  d <- ncol(Sigma.u)
  cnt <- matrix(0, n, n.mc)
  U.hist <- matrix(0, n, d*n.mc)
  X.tilde <- cbind(1, X)
  mean.u <- X.tilde%*%B.tilde
  Chol.u <- t(chol(Sigma.u))
  
  # Draw U ----
  for (idx in 1:n.mc) {
    U <- matrix(0, n, d)
    for (i in 1:n) {
      cntt <- 0
      while (TRUE) {
        cntt <- cntt + 1
        if (!uniform) {
          U.prime <- c(mean.u[i, ] + Chol.u%*%rnorm(d))
        } else {
          U.prime <- runif(d, u.lb, u.ub)
        }
        c.prime <- f(U.prime[1], U.prime[2])$c
        if (c.prime == c[i]) {
          U[i, ] <- U.prime
          break
        }
      }
      cnt[i, idx] <- cntt
    }
    U.hist[, ((idx-1)*d+1):(idx*d)] <- U
  }
  
  return(list(U.hist = U.hist, cnt = cnt))
}

#------------- STICK-BREAKING ---------------
BinaryEncode <- function (c, m) {
  # Binary encoding of categorical variables
  n <- length(c)
  n.mat <- matrix(0, n, m)
  N.mat <- matrix(0, n, m) 
  for (i in 1:n) {
    n.mat[i, c[i]] <- 1
    N.mat[i, 1:c[i]] <- 1
  }
  n.mat <- n.mat[, -m]
  N.mat <- N.mat[, -m]
  kappa.mat <- n.mat - N.mat/2
  return(list(n.mat = n.mat, N.mat = N.mat, kappa.mat = kappa.mat))
}

StickBreak <- function (prob.tilde) {
  # The mapping from (m-1) stick-breaking probabilities to m multinomial probabilities (adding up to 1).
  m <- length(prob.tilde) + 1
  prob.tilde <- c(prob.tilde, 1)
  prob <- rep(0, m)
  prob[1] <- prob.tilde[1]
  for (j in 2:m) {
    prob[j] <- prod(1 - prob.tilde[1:(j-1)])*prob.tilde[j]
  }
  return(prob)
}

StickBreak.inv <- function (prob) {
  # The mapping from m multinomial probabilities (adding up to 1) to (m-1) stick-breaking probabilities
  m <- length(prob)
  prob.tilde <- rep(0, m-1)
  prob.tilde[1] <- prob[1]
  for (j in 2:(m-1)){
    prob.tilde[j] <- prob[j]/(1 - sum(prob[1:(j-1)]))
  }
  return(prob.tilde)
}

Classify.stickbreak <- function (gamma.mat, U) {
  # Calculate the stick-breaking and multinomial probabilities for given U from gamma
  n <- nrow(U)
  d <- ncol(U)
  m <- nrow(gamma.mat) + 1
  ##
  U.tilde <- cbind(1, U)
  psi.mat <- U.tilde%*%t(gamma.mat)
  prob.tilde <- 1/(1 + exp(-psi.mat))
  prob <- matrix(0, n, m)
  for (i in 1:n) {
    prob[i, ] <- StickBreak(prob.tilde[i, ]) 
  }
  ##
  return(list(psi.mat = psi.mat, prob.tilde = prob.tilde, prob = prob))
}

#------------- SIMULATION ---------------
Branin <- function (x1, x2, 
                    x1.ub = 2, x1.lb = -2, x2.ub = 2, x2.lb = -2, 
                    a = 1, b = 5.1/(4*pi^2), c = 5/pi, s = 10, r = 6, t = 1/(8*pi)) {
  # rescale 
  x1 <- matrix((x1 - x1.lb)/(x1.ub - x1.lb)*(10 - (-5)) + (-5), ncol = 1)  #[-5, 10]
  x2 <- matrix((x2 - x2.lb)/(x2.ub - x2.lb)*(15 - 0) + 0, ncol = 1)  #[0, 15]
  x <- cbind(x1, x2)
  
  # f value
  f <- c(a*(x2 - b*x1^2 + c*x1 - r)^2 + s*(1-t)*cos(x1) + s)
  
  # linear boundary
  # line 1
  p11 <- c(10, 8)
  p12 <- c(-1, 15)
  slope1 <- (p12[2] - p11[2])/(p12[1] - p11[1])
  intp1 <- p11[2] - slope1*p11[1]
  # line 2
  p21 <- c(-5, 10)
  p22 <- c(0, 0)
  slope2 <- (p22[2] - p21[2])/(p22[1] - p21[1])
  intp2 <- p21[2] - slope2*p21[1]

  # classify
  d1 <- x2 - (intp1 + slope1*x1)
  d2 <- x2 - (intp2 + slope2*x1)
  c <- rep(0, length(f))
  c[d1 > 0] <- 1
  c[(d1 <= 0) & (d2 < 0)] <- 2
  
  # Voronoi diagram
  x.min1 <- c(-pi, 12.275)
  x.min2 <- c(pi, 2.275)
  x.min3 <- c(9.42478, 2.475)
  d11 <- apply(x, 1, function(x) sum((x-x.min1)^2))
  d22 <- apply(x, 1, function(x) sum((x-x.min2)^2))
  d33 <- apply(x, 1, function(x) sum((x-x.min3)^2))
  
  c[(d1 <= 0) & (d2 >= 0)] <- 2 + apply(cbind(d11,d22,d33)[(d1 <= 0) & (d2 >= 0), ,drop = FALSE], 1, which.min)
  
  # return
  return(list(f = f, c = c))
}

xexp <- function (x1, x2) {
  x <- cbind(x1, x2)
  
  f <- x1*exp(-x1^2 - x2^2)
  c <- rep(0, length(f))
  c[x2 >= 1.3] <- 1
  c[(x2 < 1.3) & (x1 <= 0)] <- 2
  c[(x2 < 1.3) & (x1 > 0)] <- 3
  
  # return
  return(list(f = f, c = c))
}

GetInitialU <- function (n, m, c, u.obs.idx, U.obs) {
# This function generates initial U according to a normal model fitted from the observed U.
    d <- ncol(U.obs)
  U <- matrix(0, n, d)
  for (j in 1:m) {
    U.obs.j <- U.obs[c[u.obs.idx] == j, ]
    mu.j <- apply(U.obs.j, 2, mean)
    U.obs.j <- t(t(U.obs.j) - mu.j)
    Sigma.j <-  t(U.obs.j)%*%U.obs.j/sum(c[u.obs.idx] == j)
    U[c == j, ] <- rmvnorm(sum(c == j), mu.j, Sigma.j)
  }
  U[u.obs.idx, ] <- U.obs
  return(U)
}



