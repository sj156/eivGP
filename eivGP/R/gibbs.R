# Generated from _main.Rmd: do not edit by hand

#' Unified Gibbs Sampler for EIV-GP Inference
#' 
#' @param input.type String. "ordinal" or "nominal".
#' 
#' @param Gibbs.param List of initial parameter values. Depending on `input.type`, it must contain:
#'   \itemize{
#'     \item beta0 : numeric vector of length p+1. Initial regression coefficients (intercept + p covariates). \strong{(ordinal only)}
#'     \item sigma.u20 : positive scalar. Initial value for latent variable variance \eqn{\sigma_u^2}. \strong{(ordinal only)}
#'     \item sigma.error20 : positive scalar. Initial value for error variance \eqn{\sigma_\epsilon^2}. \strong{(both)}
#'     \item Tau0 : numeric vector of length J-1. Initial thresholds (increasing order) for J categories. \strong{(ordinal only)}
#'     \item U0 : Initial latent variables. In ordinal: numeric vector of length n; in nominal: n x d matrix. \strong{(both)}
#'     \item B.tilde0 : numeric matrix of dimension (p+1) x d. Initial coefficient matrix \eqn{\tilde{B}}. \strong{(nominal only)}
#'     \item Sigma.u0 : positive definite d x d matrix. Initial covariance matrix \eqn{\Sigma_u}. \strong{(nominal only)}
#'     \item gamma.mat0 : numeric matrix of dimension (m-1) x (d+1). Initial threshold matrix \eqn{\Gamma}. \strong{(nominal only)}
#'   }
#' 
#' @param Gibbs.data List of data inputs. Common and model-specific elements:
#'   \itemize{
#'     \item X : numeric matrix (n x p). Covariates (without intercept). \strong{(both)}
#'     \item c : integer vector of length n. Observed category labels (1 to J, or 1 to m). \strong{(both)}
#'     \item y : numeric vector of length n. Continuous response. \strong{(both)}
#'     \item u.obs.idx : integer vector. Indices of observations with known U values (can be empty). \strong{(both)}
#'     \item U.obs : Observed latent values for u.obs.idx. Vector in ordinal, matrix in nominal. \strong{(both)}
#'     \item u.lb : numeric scalar or -Inf. Lower bound for U. \strong{(ordinal only)}
#'     \item u.ub : numeric scalar or Inf. Upper bound for U. \strong{(ordinal only)}
#'     \item m : integer. Number of categories in nominal response. \strong{(nominal only)}
#'   }
#' 
#' @param Gibbs.hyper List of hyperparameters. Model-specific elements:
#'   \itemize{
#'     \item v0 : prior variance for intercept (beta). \strong{(ordinal only)}
#'     \item g : scaling factor for prior precision on coefficients (beta). \strong{(ordinal only)}
#'     \item IG.a.u : shape for Inverse-Gamma prior on \eqn{\sigma_u^2}. \strong{(ordinal only)}
#'     \item IG.b.u : rate for Inverse-Gamma prior on \eqn{\sigma_u^2}. \strong{(ordinal only)}
#'     \item IG.a.error : shape for Inverse-Gamma prior on \eqn{\sigma_\epsilon^2}. \strong{(both)}
#'     \item IG.b.error : rate for Inverse-Gamma prior on \eqn{\sigma_\epsilon^2}. \strong{(both)}
#'     \item nu0 : degrees of freedom for Inverse-Wishart prior on \eqn{\Sigma_u}. \strong{(nominal only)}
#'     \item Psi0 : positive definite scale matrix for Inverse-Wishart prior on \eqn{\Sigma_u}. \strong{(nominal only)}
#'     \item q0 : prior parameter related to Polya-Gamma augmentation. \strong{(nominal only)}
#'     \item q.gamma : prior parameter related to Polya-Gamma augmentation. \strong{(nominal only)}
#'   }
#' 
#' @param GP.param List of GP parameters.
#'   \itemize{
#'     \item theta : numeric vector. Correlation length parameters (one per covariate + one for latent U). \strong{(both)}
#'     \item rho : positive scalar. Signal-to-noise ratio: \eqn{\rho = \sigma_f / \sigma_\epsilon}. \strong{(both)}
#'   }
#' 
#' @param Gibbs.other List of other MCMC settings.
#'   \itemize{
#'     \item num.iter : integer. Number of MCMC iterations. \strong{(both)}
#'     \item lmd : positive scalar. Step size for the MALA update of U. \strong{(both)}
#'     \item simu.f : logical. If TRUE, sample the latent GP function f at each iteration. \strong{(both)}
#'     \item design.map : logical. (Reserved). \strong{(both)}
#'     \item seed : integer or NULL. Random seed for reproducibility. \strong{(both)}
#'   }
#' 
#' @return A list of posterior samples (hist = history matrices/vectors).
#' \describe{
#'   \item{For ordinal model:}{
#'     \itemize{
#'       \item beta.hist : (p+1) x num.iter matrix. Each column is a draw of the regression coefficients.
#'       \item sigma.u2.hist : numeric vector of length num.iter. Draws of \eqn{\sigma_u^2}.
#'       \item sigma.error2.hist : numeric vector of length num.iter. Draws of \eqn{\sigma_{\epsilon}^2}.
#'       \item Tau.hist : (J-1) x num.iter matrix. Each column is a draw of the thresholds \eqn{\tau}.
#'       \item U.hist : n x num.iter matrix. Each column is a draw of the latent variables \eqn{U}.
#'       \item accep : logical vector of length num.iter. Acceptance indicators for the MALA update of U.
#'       \item f.hist : n x num.iter matrix (only when simu.f = TRUE). Draws of the latent GP function \eqn{f}.
#'     }
#'   }
#'   \item{For nominal model:}{
#'     \itemize{
#'       \item B.tilde.hist : (p+1) x (d * num.iter) matrix. Each block of d consecutive columns is a draw of \eqn{\tilde{B}}.
#'       \item Sigma.u.hist : d x (d * num.iter) matrix. Each block of d columns is the d x d \eqn{\Sigma_u} stored column-wise.
#'       \item sigma.error2.hist : numeric vector of length num.iter. Draws of \eqn{\sigma_{\epsilon}^2}.
#'       \item gamma.mat.hist : (m-1) x ((d+1)*num.iter) matrix. Each block of (d+1) columns is a draw of \eqn{\Gamma}.
#'       \item U.hist : n x (d * num.iter) matrix. Each block of d columns is a draw of the latent matrix \eqn{U} (n x d).
#'       \item accep : logical vector of length num.iter. Acceptance indicators.
#'       \item f.hist : n x num.iter matrix (if simu.f = TRUE). Draws of the latent function \eqn{f}.
#'     }
#'   }
#' }
#' @export
Gibbs <- function(input.type, Gibbs.param, Gibbs.data, Gibbs.hyper, GP.param, Gibbs.other) {
  
  if (!(input.type %in% c("ordinal", "nominal"))) {
    stop("Input Error: 'input.type' must be 'ordinal' or 'nominal'.")
  }
  
  # Set random seed if provided in Gibbs.other
  if (!is.null(Gibbs.other$seed)) {
    set.seed(Gibbs.other$seed)
  }
  
  if (input.type == "ordinal") {
    message("Running Gibbs Sampler for Ordinal model...")
    return(Gibbs.ordinal(Gibbs.param, Gibbs.data, Gibbs.hyper, GP.param, Gibbs.other))
  } else {
    message("Running Gibbs Sampler for Nominal model...")
    return(Gibbs.nominal(Gibbs.param, Gibbs.data, Gibbs.hyper, GP.param, Gibbs.other))
  }
}

#' Gibbs sampler for ordinal case
#' @param Gibbs.param: see above
#' @param Gibbs.data: see above 
#' @param Gibbs.hyper: see above
#' @param GP.param: see above
#' @param Gibbs.other: see above
#' @return see above
Gibbs.ordinal <- function (Gibbs.param, Gibbs.data, Gibbs.hyper, GP.param, Gibbs.other) {
  
  # 1. Unpack Inputs
  
  # Gibbs.other
  num.iter <- Gibbs.other$num.iter
  lmd <- Gibbs.other$lmd
  simu.f <- if(!is.null(Gibbs.other$simu.f)) Gibbs.other$simu.f else TRUE
  design.map <- if(!is.null(Gibbs.other$design.map)) Gibbs.other$design.map else FALSE
  
  # Gibbs.param
  beta0 <- Gibbs.param$beta0
  sigma.u20 <- Gibbs.param$sigma.u20
  sigma.error20 <- Gibbs.param$sigma.error20
  Tau0 <- Gibbs.param$Tau0
  U0 <- Gibbs.param$U0
  
  # Gibbs.data
  X <- Gibbs.data$X; c <- Gibbs.data$c; y <- Gibbs.data$y
  u.obs.idx <- Gibbs.data$u.obs.idx; U.obs <- Gibbs.data$U.obs
  u.lb <- if(!is.null(Gibbs.data$u.lb)) Gibbs.data$u.lb else -Inf
  u.ub <- if(!is.null(Gibbs.data$u.ub)) Gibbs.data$u.ub else Inf
  
  # Gibbs.hyper
  v0 <- Gibbs.hyper$v0; g <- Gibbs.hyper$g
  IG.a.u <- Gibbs.hyper$IG.a.u; IG.b.u <- Gibbs.hyper$IG.b.u
  IG.a.error <- Gibbs.hyper$IG.a.error; IG.b.error <- Gibbs.hyper$IG.b.error
  
  # GP.param
  theta <- GP.param$theta; rho <- GP.param$rho
  
  # 2. Compute Constants and Initialize
  n <- length(c); m.categories <- length(Tau0) + 1; d.x <- ncol(X); n.obs <- length(u.obs.idx)
  u.mis.idx <- if (n.obs > 0) c(1:n)[-u.obs.idx] else c(1:n)
  
  X.tilde <- cbind(1, X); XX.tilde <- t(X.tilde) %*% X.tilde
  V0 <- matrix(0, d.x+1, d.x+1); V0[1,1] <- v0; V0[-1,-1] <- g * XX.tilde[-1,-1]
  Vn <- V0 + XX.tilde; Vn.inv <- solve(Vn)
  
  # Allocate storage
  beta.hist <- matrix(0, d.x + 1, num.iter); sigma.u2.hist <- rep(0, num.iter)
  sigma.error2.hist <- rep(0, num.iter); Tau.hist <- matrix(0, m.categories-1, num.iter)
  U.hist <- matrix(0, n, num.iter); accep <- rep(0, num.iter); f.hist <- matrix(0, n, num.iter)
  
  beta <- beta0; sigma.u2 <- sigma.u20; sigma.error2 <- sigma.error20; Tau <- Tau0; U <- U0
  if (n.obs > 0) U[u.obs.idx] <- U.obs
  
  Theta <- c(theta, rho^2 * sigma.error2)
  #print(ncol(cbind(X,U))+1)
  #print(length(Theta))
  KWW <- KernelMatrix(cbind(X, U), Theta); C <- KWW + sigma.error2 * diag(n)
  C.chol <- t(chol(C)); r <- backsolve(t(C.chol), forwardsolve(C.chol, y))
  
  # 3. MCMC Loop
  for (k in 1:num.iter) {
    # Update beta and sigma.u2
    mu.n <- Vn.inv %*% t(X.tilde) %*% U
    sigma.u2 <- 1/rgamma(1, shape = IG.a.u + n/2, rate = IG.b.u + (t(U) %*% U - t(mu.n) %*% Vn %*% mu.n)/2)
    beta <- c(rmvnorm(1, mu.n, sigma.u2 * Vn.inv))
    sigma.u2.hist[k] <- sigma.u2; beta.hist[, k] <- beta 
    
    # Update sigma.error2
    sigma.error2.prime <- 1/rgamma(1, shape = IG.a.error + n/2, rate = IG.b.error + sigma.error2/2 * t(y) %*% r)
    Theta <- c(theta, rho^2 * sigma.error2.prime)
    KWW <- (sigma.error2.prime/sigma.error2) * KWW; C <- (sigma.error2.prime/sigma.error2) * C
    C.chol <- sqrt(sigma.error2.prime/sigma.error2) * C.chol; r <- (sigma.error2/sigma.error2.prime) * r
    sigma.error2 <- sigma.error2.prime; sigma.error2.hist[k] <- sigma.error2
    
    # Update Tau
    for (j in 1:(m.categories-1)) Tau[j] <- runif(1, max(U[c == j]), min(U[c == j+1]))
    Tau.hist[, k] <- Tau
    
    # Update U via LMC
    phi <- backsolve(t(C.chol), rnorm(n))
    GU <- Gradient.ordinal(U, r, phi, KWW, theta[d.x + 1], X.tilde, beta, sigma.u2)
    
    Tau.aug <- c(u.lb, Tau, u.ub); U.mis <- U[u.mis.idx]
    U.mis.prime <- DrawTrunNormal(U.mis + lmd^2/2 * GU[u.mis.idx], rep(lmd, length(u.mis.idx)), 
                                  Tau.aug[c[u.mis.idx]], Tau.aug[c[u.mis.idx]+1])
    
    U.prime <- U; U.prime[u.mis.idx] <- U.mis.prime
    KWW.prime <- KernelMatrix(cbind(X, U.prime), Theta); C.prime <- KWW.prime + sigma.error2 * diag(n)
    C.prime.chol <- t(chol(C.prime)); r.prime <- backsolve(t(C.prime.chol), forwardsolve(C.prime.chol, y))
    
    GU.prime <- Gradient.ordinal(U.prime, r.prime, phi, KWW.prime, theta[d.x + 1], X.tilde, beta, sigma.u2)
    
    # Accept-reject step
    if(log(runif(1)) <= 0) {  
      U <- U.prime; KWW <- KWW.prime; C <- C.prime; C.chol <- C.prime.chol; r <- r.prime; accep[k] <- 1
    }
    
    U.hist[, k] <- U
    if (simu.f) {
      LiK <- apply(KWW, 2, function(b) forwardsolve(C.chol, b))
      f.hist[, k] <- rmvnorm(1, KWW %*% r, KWW - t(LiK) %*% LiK)
    }
  }
  return(list(beta.hist = beta.hist, sigma.u2.hist = sigma.u2.hist, sigma.error2.hist = sigma.error2.hist,
              Tau.hist = Tau.hist, U.hist = U.hist, accep = accep, f.hist = f.hist))
}


#' Gibbs sampler for nominal case (Note that well organized method now)
#' @param Gibbs.param: see above
#' @param Gibbs.data: see above 
#' @param Gibbs.hyper: see above
#' @param GP.param: see above
#' @param Gibbs.other: see above
#' @return see above
Gibbs.nominal <- function (Gibbs.param, Gibbs.data, Gibbs.hyper, GP.param, Gibbs.other) {
  
  # 1. Unpack Inputs
  
  # Gibbs.other
  num.iter <- Gibbs.other$num.iter
  lmd <- Gibbs.other$lmd
  simu.f <- if(!is.null(Gibbs.other$simu.f)) Gibbs.other$simu.f else TRUE
  design.map <- if(!is.null(Gibbs.other$design.map)) Gibbs.other$design.map else FALSE
  
  # Gibbs.param
  B.tilde0 <- Gibbs.param$B.tilde0
  Sigma.u0 <- Gibbs.param$Sigma.u0
  sigma.error20 <- Gibbs.param$sigma.error20
  gamma.mat0 <- Gibbs.param$gamma.mat0
  U0 <- Gibbs.param$U0
  
  # Gibbs.data
  X <- Gibbs.data$X; c <- Gibbs.data$c; y <- Gibbs.data$y; m <- Gibbs.data$m
  u.obs.idx <- Gibbs.data$u.obs.idx; U.obs <- Gibbs.data$U.obs
  
  # Gibbs.hyper
  v0 <- Gibbs.hyper$v0; g <- Gibbs.hyper$g
  nu0 <- Gibbs.hyper$nu0; Psi0 <- Gibbs.hyper$Psi0
  q0 <- Gibbs.hyper$q0; q.gamma <- Gibbs.hyper$q.gamma
  IG.a.error <- Gibbs.hyper$IG.a.error; IG.b.error <- Gibbs.hyper$IG.b.error
  
  # GP.param
  theta <- GP.param$theta; rho <- GP.param$rho
  
  # 2. Compute Constants and Initialize
  n <- nrow(X); d.x <- ncol(X); d <- ncol(U0); n.obs <- length(u.obs.idx)
  u.mis.idx <- if (n.obs > 0) c(1:n)[-u.obs.idx] else c(1:n)
  
  bin.encoding <- BinaryEncode(c, m); N.mat <- bin.encoding$N.mat; kappa.mat <- bin.encoding$kappa.mat
  X.tilde <- cbind(1, X); XX.tilde <- t(X.tilde) %*% X.tilde
  Vn.inv <- as.symmetric.matrix(solve(matrix(0, d.x+1, d.x+1) + XX.tilde))
  
  # Allocate storage
  B.tilde.hist <- matrix(0, d.x+1, d*num.iter); Sigma.u.hist <- matrix(0, d, d*num.iter)
  sigma.error2.hist <- rep(0, num.iter); gamma.mat.hist <- matrix(0, m-1, (d+1)*num.iter)
  U.hist <- matrix(0, n, d*num.iter); accep <- rep(0, num.iter); f.hist <- matrix(0, n, num.iter)
  
  B.tilde <- B.tilde0; Sigma.u <- Sigma.u0; sigma.error2 <- sigma.error20; gamma.mat <- gamma.mat0; U <- U0
  if (n.obs > 0) U[u.obs.idx, ] <- U.obs
  omega.mat <- matrix(0, n, m-1); mn.u <- matrix(0, n, d)
  Qn.u <- array(0, c(d,d,n))
  
  Theta <- c(theta, rho^2 * sigma.error2)
  KWW <- KernelMatrix(cbind(X, U), Theta); C <- KWW + sigma.error2 * diag(n)
  C.chol <- t(chol(C)); r <- backsolve(t(C.chol), forwardsolve(C.chol, y))
  
  # 3. MCMC Loop
  for (k.iter in 1:num.iter) {
    # Update B.tilde and Sigma.u
    Mn <- Vn.inv %*% t(X.tilde) %*% U; Psin <- Psi0 + t(U) %*% U - t(Mn) %*% Vn.inv %*% Mn 
    Sigma.u <- rinvwishart(nu0 + n, Psin); B.tilde <- rmatrixnorm(Mn, Vn.inv, Sigma.u)
    Sigma.u.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- Sigma.u
    B.tilde.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- B.tilde
    
    # Update sigma.error2
    sigma.error2.prime <- 1/rgamma(1, shape = IG.a.error + n/2, rate = IG.b.error + sigma.error2/2 * t(y) %*% r)
    Theta <- c(theta, rho^2 * sigma.error2.prime)
    KWW <- (sigma.error2.prime/sigma.error2)*KWW; C <- (sigma.error2.prime/sigma.error2)*C
    C.chol <- sqrt(sigma.error2.prime/sigma.error2)*C.chol; r <- (sigma.error2/sigma.error2.prime)*r
    sigma.error2 <- sigma.error2.prime; sigma.error2.hist[k.iter] <- sigma.error2
    
    # Update auxiliary variables (Polya-Gamma, omega, etc. logic remains unchanged)...
    
    phi <- backsolve(t(C.chol), rnorm(n))
    Sigma.u.inv <- solve(Sigma.u)
    for (i in 1:n) {
      Qn.u[,,i] <-  Sigma.u.inv + t(gamma.mat[,-1]) %*% diag(omega.mat[i,]) %*% gamma.mat[,-1]
      mn.u[i, ] <- solve(Qn.u[,,i]) %*% (Sigma.u.inv %*% t(B.tilde) %*% X.tilde[i,] + t(gamma.mat[,-1]) %*% (kappa.mat[i,] - omega.mat[i,]*gamma.mat[,1]))
    }
    
    GU <- Gradient.nominal(U, phi, r, KWW, mn.u, Qn.u, theta[d.x+1])
    U.mis <- U[u.mis.idx, ]; U.mis.prime <- U.mis + lmd^2/2 * GU[u.mis.idx, ] + lmd * rnorm(length(u.mis.idx)*d)
    
    U.prime <- U; U.prime[u.mis.idx, ] <- U.mis.prime
    KWW.prime <- KernelMatrix(cbind(X, U.prime), Theta); C.prime <- KWW.prime + sigma.error2 * diag(n)
    C.prime.chol <- t(chol(C.prime)); r.prime <- backsolve(t(C.prime.chol), forwardsolve(C.prime.chol, y))
    
    GU.prime <- Gradient.nominal(U.prime, phi, r.prime, KWW.prime, mn.u, Qn.u, theta[d.x+1])
    
    # Accept-reject step
    if (log(runif(1)) <= 0) { 
      U <- U.prime; KWW <- KWW.prime; C <- C.prime; C.chol <- C.prime.chol; r <- r.prime; accep[k.iter] <- 1
    }
    U.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- U
    
    if (simu.f) {
      LiK <- apply(KWW, 2, function(b) forwardsolve(C.chol, b))
      f.hist[, k.iter] <- rmvnorm(1, KWW %*% r, KWW - t(LiK) %*% LiK)
    }
  }
  return(list(B.tilde.hist = B.tilde.hist, Sigma.u.hist = Sigma.u.hist, sigma.error2.hist = sigma.error2.hist,
              gamma.mat.hist = gamma.mat.hist, U.hist = U.hist, accep = accep, f.hist = f.hist))
}


#' Unified Empirical Bayes Gibbs Sampler
#' 
#' @description
#' Runs a Gibbs sampler for ordinal or nominal probit models without a Gaussian process.
#' This version is designed for hyperparameter tuning within an empirical Bayes framework,
#' omitting the continuous response \eqn{y} and the associated error variance \eqn{\sigma_\epsilon^2}.
#' 
#' @param input.type String. `"ordinal"` or `"nominal"`.
#' 
#' @param Gibbs.param List of initial parameter values. Elements depend on `input.type`:
#'   \itemize{
#'     \item beta0 : numeric vector of length p+1. Initial regression coefficients (intercept + p covariates). \strong{(ordinal only)}
#'     \item sigma.u20 : positive scalar. Initial value for latent variable variance \eqn{\sigma_u^2}. \strong{(ordinal only)}
#'     \item Tau0 : numeric vector of length J-1. Initial thresholds (increasing) for J categories. \strong{(ordinal only)}
#'     \item U0 : Initial latent variables. \strong{(both)} \cr
#'           In ordinal: numeric vector of length n. \cr
#'           In nominal: numeric matrix of dimension n \eqn{\times} d.
#'     \item B.tilde0 : numeric matrix of dimension (p+1) \eqn{\times} d. Initial coefficient matrix \eqn{\tilde{B}}. \strong{(nominal only)}
#'     \item Sigma.u0 : positive definite d \eqn{\times} d matrix. Initial covariance matrix \eqn{\Sigma_u}. \strong{(nominal only)}
#'     \item gamma.mat0 : numeric matrix of dimension (m-1) \eqn{\times} (d+1). Initial threshold matrix \eqn{\Gamma}. \strong{(nominal only)}
#'   }
#' 
#' @param Gibbs.data List of data inputs. Elements depend on `input.type`:
#'   \itemize{
#'     \item X : numeric matrix (n \eqn{\times} p). Covariates (without intercept). \strong{(both)}
#'     \item c : integer vector of length n. Observed category labels (1 to J or 1 to m). \strong{(both)}
#'     \item u.obs.idx : integer vector. Indices of observations with known latent U values (can be empty). \strong{(both)}
#'     \item U.obs : Observed latent values for `u.obs.idx`. \strong{(both)} \cr
#'           In ordinal: numeric vector. \cr
#'           In nominal: numeric matrix with rows corresponding to observations.
#'     \item u.lb : numeric scalar or `-Inf`. Lower bound for the latent variable U. \strong{(ordinal only)}
#'     \item u.ub : numeric scalar or `Inf`. Upper bound for the latent variable U. \strong{(ordinal only)}
#'     \item m : integer. Number of categories in the nominal response. \strong{(nominal only)}
#'   }
#' 
#' @param Gibbs.hyper List of hyperparameters. Elements depend on `input.type`:
#'   \itemize{
#'     \item v0 : prior variance for the intercept (used in the prior covariance of \eqn{\beta} or \eqn{\tilde{B}}). \strong{(both)}
#'     \item g : scaling factor for the prior precision matrix of the regression coefficients. \strong{(both)}
#'     \item IG.a.u : shape for Inverse-Gamma prior on \eqn{\sigma_u^2}. \strong{(ordinal only)}
#'     \item IG.b.u : rate for Inverse-Gamma prior on \eqn{\sigma_u^2}. \strong{(ordinal only)}
#'     \item nu0 : degrees of freedom for the Inverse-Wishart prior on \eqn{\Sigma_u}. \strong{(nominal only)}
#'     \item Psi0 : positive definite scale matrix for the Inverse-Wishart prior on \eqn{\Sigma_u}. \strong{(nominal only)}
#'     \item q0 : prior precision for the intercept column of \eqn{\Gamma}. \strong{(nominal only)}
#'     \item q.gamma : prior precision for the slope columns of \eqn{\Gamma}. \strong{(nominal only)}
#'   }
#' 
#' @param Gibbs.other List of other MCMC settings.
#'   \itemize{
#'     \item n.mc : integer. Number of MCMC iterations. \strong{(both)}
#'     \item seed : integer or `NULL`. Random seed for reproducibility (optional). \strong{(both)}
#'   }
#' 
#' @return A list of posterior samples (history matrices/vectors).
#' \describe{
#'   \item{For ordinal model:}{
#'     \itemize{
#'       \item beta.hist : (p+1) \eqn{\times} n.mc matrix. Each column is a draw of the regression coefficients.
#'       \item sigma.u2.hist : numeric vector of length n.mc. Draws of \eqn{\sigma_u^2}.
#'       \item Tau.hist : (J-1) \eqn{\times} n.mc matrix. Each column is a draw of the thresholds \eqn{\tau}.
#'       \item U.hist : n \eqn{\times} n.mc matrix. Each column is a draw of the latent variables U.
#'     }
#'   }
#'   \item{For nominal model:}{
#'     \itemize{
#'       \item B.tilde.hist : (p+1) \eqn{\times} (d * n.mc) matrix. Each block of d consecutive columns is a draw of \eqn{\tilde{B}}.
#'       \item Sigma.u.hist : d \eqn{\times} (d * n.mc) matrix. Each block of d columns is the d \eqn{\times} d \eqn{\Sigma_u} stored column‑wise.
#'       \item gamma.mat.hist : (m-1) \eqn{\times} ((d+1) * n.mc) matrix. Each block of (d+1) columns is a draw of \eqn{\Gamma}.
#'       \item U.hist : n \eqn{\times} (d * n.mc) matrix. Each block of d columns is a draw of the latent matrix U (n \eqn{\times} d).
#'     }
#'   }
#' }
#' 
#' @export
Gibbs.EB <- function(input.type, Gibbs.param, Gibbs.data, Gibbs.hyper, Gibbs.other) {
  
  if (!(input.type %in% c("ordinal", "nominal"))) {
    stop("Input Error: 'input.type' must be 'ordinal' or 'nominal'.")
  }
  
  # Set random seed if provided in Gibbs.other
  if (!is.null(Gibbs.other$seed)) {
    set.seed(Gibbs.other$seed)
  }
  
  if (input.type == "ordinal") {
    message("Running Empirical Bayes Sampler for Ordinal model...")
    return(Gibbs.EB.ordinal(Gibbs.param, Gibbs.data, Gibbs.hyper, Gibbs.other))
  } else {
    message("Running Empirical Bayes Sampler for Nominal model...")
    return(Gibbs.EB.nominal(Gibbs.param, Gibbs.data, Gibbs.hyper, Gibbs.other))
  }
}


#' Empirical Bayes for Gibbs for tuning, Ordinal case
#' 
#' @param Gibbs.param see above
#' @param Gibbs.data see above
#' @param Gibbs.hyper see above
#' @param Gibbs.other see above
#' @return related list (for detail, see above)
Gibbs.EB.ordinal <- function (Gibbs.param, Gibbs.data, Gibbs.hyper, Gibbs.other) {
  
  # 1. Unpack Inputs
  
  # Gibbs.other
  n.mc <- Gibbs.other$n.mc
  
  # Gibbs.param
  beta0 <- Gibbs.param$beta0
  sigma.u20 <- Gibbs.param$sigma.u20
  Tau0 <- Gibbs.param$Tau0
  U0 <- Gibbs.param$U0
  
  # Gibbs.data
  X <- Gibbs.data$X 
  c <- Gibbs.data$c
  u.obs.idx <- Gibbs.data$u.obs.idx
  U.obs <- Gibbs.data$U.obs
  
  # Gibbs.hyper
  v0 <- Gibbs.hyper$v0
  g <- Gibbs.hyper$g
  IG.a.u <- Gibbs.hyper$IG.a.u
  IG.b.u <- Gibbs.hyper$IG.b.u
  
  # 2. Compute Constants and Initialize
  n <- nrow(X); m.categories <- length(Tau0) + 1; d.x <- ncol(X)
  X.tilde <- cbind(1, X); XX.tilde <- t(X.tilde) %*% X.tilde
  
  V0 <- matrix(0, d.x+1, d.x+1); V0[1,1] <- v0; V0[-1,-1] <- g * XX.tilde[-1,-1]
  Vn <- V0 + XX.tilde; Vn.inv <- solve(Vn)
  
  u.mis.idx <- if (length(u.obs.idx) > 0) c(1:n)[-u.obs.idx] else c(1:n)
  
  beta <- beta0; sigma.u2 <- sigma.u20; Tau <- Tau0; U <- U0
  if (length(u.obs.idx) > 0) U[u.obs.idx] <- U.obs 
  
  # Allocate storage
  beta.hist <- matrix(0, d.x+1, n.mc); sigma.u2.hist <- rep(0, n.mc)
  Tau.hist <- matrix(0, m.categories-1, n.mc); U.hist <- matrix(0, n, n.mc)
  
  # 3. MCMC Loop
  for (k in 1:n.mc) {
    # Update beta & sigma.u2
    mu.n <- Vn.inv %*% t(X.tilde) %*% U
    sigma.u2 <- 1/rgamma(1, shape = IG.a.u + n/2, rate = IG.b.u + (t(U) %*% U - t(mu.n) %*% Vn %*% mu.n)/2)
    beta <- c(rmvnorm(1, mu.n, sigma.u2 * Vn.inv))
    
    sigma.u2.hist[k] <- sigma.u2; beta.hist[, k] <- beta 
    
    # Update Tau
    for (j in 1:(m.categories-1)) {
      Tau[j] <- runif(1, max(U[c == j]), min(U[c == j+1]))
    }
    Tau.hist[, k] <- Tau
    
    # Update missing U
    Tau.aug <- c(u.lb, Tau, u.ub)
    U[u.mis.idx] <- DrawTrunNormal(X.tilde[u.mis.idx, ] %*% beta, 
                                   rep(sqrt(sigma.u2), length(u.mis.idx)), 
                                   Tau.aug[c[u.mis.idx]], Tau.aug[c[u.mis.idx]+1])
    U.hist[, k] <- U
  }
  
  return(list(beta.hist = beta.hist, 
              sigma.u2.hist = sigma.u2.hist,
              Tau.hist = Tau.hist,
              U.hist = U.hist))
}


#' Empirical Bayes for Gibbs for tuning, Nominal case
#' 
#' @param Gibbs.param see above
#' @param Gibbs.data see above
#' @param Gibbs.hyper see above
#' @param Gibbs.other see above
#' @return related list (for detail, see above)
Gibbs.EB.nominal <- function (Gibbs.param, Gibbs.data, Gibbs.hyper, Gibbs.other) {
  
  # 1. Unpack Inputs
  
  # Gibbs.other
  n.mc <- Gibbs.other$n.mc
  
  # Gibbs.param
  B.tilde0 <- Gibbs.param$B.tilde0
  Sigma.u0 <- Gibbs.param$Sigma.u0
  gamma.mat0 <- Gibbs.param$gamma.mat0
  U0 <- Gibbs.param$U0
  
  # Gibbs.data
  X <- Gibbs.data$X
  c <- Gibbs.data$c
  m.categories <- Gibbs.data$m
  u.obs.idx <- Gibbs.data$u.obs.idx
  U.obs <- Gibbs.data$U.obs
  
  # Gibbs.hyper
  v0 <- Gibbs.hyper$v0
  g <- Gibbs.hyper$g
  nu0 <- Gibbs.hyper$nu0
  Psi0 <- Gibbs.hyper$Psi0
  q0 <- Gibbs.hyper$q0
  q.gamma <- Gibbs.hyper$q.gamma
  
  # 2. Compute Constants and Initialize
  n <- nrow(X); d.x <- ncol(X); d <- ncol(U0)
  u.mis.idx <- if (length(u.obs.idx) > 0) c(1:n)[-u.obs.idx] else c(1:n)
  
  bin.encoding <- BinaryEncode(c, m.categories)
  N.mat <- bin.encoding$N.mat; kappa.mat <- bin.encoding$kappa.mat
  
  X.tilde <- cbind(1, X); XX.tilde <- t(X.tilde) %*% X.tilde
  V0 <- matrix(0, d.x+1, d.x+1); V0[1,1] <- v0; V0[-1,-1] <- g * XX.tilde[-1,-1]
  Vn <- V0 + XX.tilde; Vn.inv <- solve(Vn)
  Q0.gamma <- diag(c(q0, rep(q.gamma, d)))
  
  B.tilde <- B.tilde0; Sigma.u <- Sigma.u0; gamma.mat <- gamma.mat0; U <- U0
  #browser()
  if (length(u.obs.idx) > 0) U[u.obs.idx, ] <- U.obs
  omega.mat <- matrix(0, n, m.categories-1)
  
  # Allocate storage
  B.tilde.hist <- matrix(0, d.x+1, d*n.mc)
  Sigma.u.hist <- matrix(0, d, d*n.mc)
  gamma.mat.hist <- matrix(0, m.categories-1, (d+1)*n.mc)
  U.hist <- matrix(0, n, d*n.mc)
  
  # 3. MCMC Loop
  for (k.iter in 1:n.mc) {
    # Update B.tilde & Sigma.u
    nun <- nu0 + n
    Mn <- Vn.inv %*% t(X.tilde) %*% U
    Psin <- Psi0 + t(U) %*% U - t(Mn) %*% Vn %*% Mn
    
    Sigma.u <- rinvwishart(nun, Psin)
    B.tilde <- rmatrixnorm(Mn, Vn.inv, Sigma.u)
    
    Sigma.u.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- Sigma.u
    B.tilde.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- B.tilde 
    
    # Update Polya-Gamma auxiliary variables
    U.tilde <- cbind(1, U)
    psi.mat <- U.tilde %*% t(gamma.mat)
    # Note: Requires 'pgdraw' package
    omega.mat[N.mat == 1] <- pgdraw(1, psi.mat[N.mat == 1])
    
    # Update nominal choice coefficients (gamma.mat)
    for (j in 1:(m.categories-1)) {
      Qn.gamma <- Q0.gamma + t(U.tilde) %*% diag(omega.mat[, j]) %*% U.tilde
      mn.gamma <- solve(Qn.gamma) %*% t(U.tilde) %*% kappa.mat[, j]
      gamma.mat[j, ] <- rmvnorm(1, mn.gamma, solve(Qn.gamma))
    }
    gamma.mat.hist[, ((k.iter-1)*(d+1)+1):(k.iter*(d+1))] <- gamma.mat
    
    # Update missing U
    Sigma.u.inv <- solve(Sigma.u)
    for (i in u.mis.idx) {
      Qn.u <- Sigma.u.inv + t(gamma.mat[,-1]) %*% diag(omega.mat[i,]) %*% gamma.mat[,-1]
      Qn.u.inv <- solve(Qn.u)
      mn.u <- Qn.u.inv %*% (Sigma.u.inv %*% t(B.tilde) %*% X.tilde[i,] + t(gamma.mat[,-1]) %*% (kappa.mat[i,] - omega.mat[i,] * gamma.mat[,1]))
      U[i, ] <- rmvnorm(1, mn.u, Qn.u.inv)
    }
    U.hist[, ((k.iter-1)*d+1):(k.iter*d)] <- U
  }
  
  return(list(B.tilde.hist = B.tilde.hist,
              Sigma.u.hist = Sigma.u.hist,
              gamma.mat.hist = gamma.mat.hist,
              U.hist = U.hist))
}


#' Unified Function to Impute New Latent Variables U
#' 
#' @description
#' Draws posterior predictive samples of the latent variable \eqn{U} for new observations,
#' given the posterior samples of model parameters from a previous MCMC run.
#' 
#' @param input.type String. `"ordinal"` or `"nominal"`.
#' 
#' @param Gibbs.data List of data for the new observations. Elements depend on `input.type`:
#'   \itemize{
#'     \item X.star : numeric matrix (n.star \eqn{\times} p). Covariates for new observations (without intercept). \strong{(ordinal only)}
#'     \item c.star : integer vector of length n.star. Observed category labels for new observations. \strong{(ordinal only)}
#'     \item u.lb : numeric scalar or `-Inf`. Lower bound for the latent variable U. \strong{(ordinal only)}
#'     \item u.ub : numeric scalar or `Inf`. Upper bound for the latent variable U. \strong{(ordinal only)}
#'     \item X : numeric matrix (n \eqn{\times} p). Original covariates (used when refining U draws). \strong{(nominal only)}
#'     \item c : integer vector of length n. Original category labels. \strong{(nominal only)}
#'     \item U0 : numeric matrix of dimension n \eqn{\times} d. Initial values for the latent U (typically the last draw from the original chain). \strong{(nominal only)}
#'   }
#' 
#' @param Gibbs.hist List containing MCMC history matrices/vectors from a fitted model. Elements depend on `input.type`:
#'   \itemize{
#'     \item beta.hist : (p+1) \eqn{\times} n.mc matrix. Posterior draws of regression coefficients. \strong{(ordinal only)}
#'     \item sigma.u2.hist : numeric vector of length n.mc. Posterior draws of \eqn{\sigma_u^2}. \strong{(ordinal only)}
#'     \item Tau.hist : (J-1) \eqn{\times} n.mc matrix. Posterior draws of thresholds. \strong{(ordinal only)}
#'     \item B.tilde.hist : (p+1) \eqn{\times} (d * n.mc) matrix. Blocked posterior draws of \eqn{\tilde{B}}. \strong{(nominal only)}
#'     \item Sigma.u.hist : d \eqn{\times} (d * n.mc) matrix. Blocked posterior draws of \eqn{\Sigma_u}. \strong{(nominal only)}
#'     \item gamma.mat.hist : (m-1) \eqn{\times} ((d+1) * n.mc) matrix. Blocked posterior draws of \eqn{\Gamma}. \strong{(nominal only)}
#'   }
#' 
#' @param Gibbs.other List of additional settings.
#'   \itemize{
#'     \item seed : integer or `NULL`. Random seed for reproducibility (optional). \strong{(both)}
#'     \item n.mc : integer. Number of posterior samples to process (usually equals the number of MCMC iterations). \strong{(nominal only)}
#'     \item num.iter : integer. Number of inner Gibbs iterations per MCMC sample (used only within nominal imputation). \strong{(nominal only)}
#'   }
#' 
#' @return
#' \itemize{
#'   \item For **ordinal**: a list containing \code{U.star.hist}, an n.star \eqn{\times} n.mc matrix of imputed latent variables.
#'   \item For **nominal**: an n \eqn{\times} (d * n.mc) matrix where each block of d columns corresponds to one imputed U draw.
#' }
#' 
#' @export
impute.u <- function(input.type, Gibbs.data, Gibbs.hist, Gibbs.other) {
  
  if (!(input.type %in% c("ordinal", "nominal"))) {
    stop("Input Error: 'input.type' must be 'ordinal' or 'nominal'.")
  }
  
  # Set random seed if provided
  if (!is.null(Gibbs.other$seed)) {
    set.seed(Gibbs.other$seed)
  }
  
  if (input.type == "ordinal") {
    message("Imputing U for Ordinal model...")
    return(impute.u.ordinal(Gibbs.data, Gibbs.hist, Gibbs.other))
  } else {
    message("Imputing U for Nominal model...")
    return(impute.u.nominal(Gibbs.data, Gibbs.hist, Gibbs.other))
  }
}


#' Impute New Latent Variables U for Ordinal Models
impute.u.ordinal <- function(Gibbs.data, Gibbs.hist, Gibbs.other) {
  
  # 1. Unpack Inputs
  
  # Gibbs.data
  X.star <- Gibbs.data$X.star
  c.star <- Gibbs.data$c.star
  u.lb <- if(!is.null(Gibbs.data$u.lb)) Gibbs.data$u.lb else -Inf
  u.ub <- if(!is.null(Gibbs.data$u.ub)) Gibbs.data$u.ub else Inf
  
  # Gibbs.hist
  beta.hist <- Gibbs.hist$beta.hist
  sigma.u2.hist <- Gibbs.hist$sigma.u2.hist
  Tau.hist <- Gibbs.hist$Tau.hist
  
  # Gibbs.other
  seed <- if(!is.null(Gibbs.other$seed)) Gibbs.other$seed else 222
  
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
  
  return(list(U.star.hist = U.star.hist))
}


#' Impute New Latent Variables U for Nominal Models
impute.u.nominal <- function(Gibbs.data, Gibbs.hist, Gibbs.other) {
  
  # 1. Unpack Inputs
  
  # Gibbs.data
  X <- Gibbs.data$X
  c <- Gibbs.data$c
  U0 <- Gibbs.data$U0
  # Gibbs.hist
  B.tilde.hist <- Gibbs.hist$B.tilde.hist
  Sigma.u.hist <- Gibbs.hist$Sigma.u.hist
  gamma.mat.hist <- Gibbs.hist$gamma.mat.hist
  
  # Gibbs.other
  n.mc <- Gibbs.other$n.mc
  seed <-  if(!is.null(Gibbs.other$seed)) Gibbs.data$seed else 222
  
  #Constant
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


#' Impute Continuous Response y via Gaussian Process Kriging
#'
#' @description
#' Predicts the continuous response \eqn{y} (and the latent GP function \eqn{f})
#' at new covariate locations, given posterior samples of the latent variables
#' \eqn{U} and error variance \eqn{\sigma_\epsilon^2}. Applies to both ordinal
#' and nominal models.
#'
#' @param Gibbs.data List of data elements.
#'   \itemize{
#'     \item y : numeric vector of length n. Original continuous response. \strong{(both)}
#'     \item X : numeric matrix (n \eqn{\times} p). Original covariates. \strong{(both)}
#'     \item X.star : numeric matrix (n.star \eqn{\times} p). New covariates at which to predict. \strong{(both)}
#'     \item d : integer. Dimension of the latent space (1 for ordinal, \eqn{d \ge 1} for nominal). \strong{(both)}
#'   }
#'
#' @param Gibbs.hist List containing MCMC history objects. Must contain:
#'   \itemize{
#'     \item sigma.error2.hist : numeric vector of length n.mc. Posterior draws of \eqn{\sigma_\epsilon^2}. \strong{(both)}
#'     \item U.hist : n \eqn{\times} (d * n.mc) matrix. Blocked posterior draws of the latent U for original observations. \strong{(both)}
#'     \item U.star.hist : Imputed latent U for new observations. \strong{(both)} \cr
#'           If \code{single.ustar = TRUE}, a fixed matrix of size n.star \eqn{\times} d; \cr
#'           otherwise, an n.star \eqn{\times} (d * n.mc) matrix, with blocks matching \code{U.hist}.
#'   }
#'
#' @param GP.param List of GP parameters.
#'   \itemize{
#'     \item theta : numeric vector. Correlation length parameters. \strong{(both)}
#'     \item rho : positive scalar. Signal-to-noise ratio \eqn{\rho = \sigma_f / \sigma_\epsilon}. \strong{(both)}
#'   }
#'
#' @param Gibbs.other List of optional settings.
#'   \itemize{
#'     \item single.ustar : logical. If \code{TRUE}, \code{U.star.hist} is treated as a single fixed U.star matrix (used for all iterations). Default \code{FALSE}. \strong{(both)}
#'     \item separate : logical. If \code{TRUE}, draws \eqn{f} and \eqn{y} independently per test point (diagonal covariance approximation). Default \code{FALSE}. \strong{(both)}
#'     \item seed : integer or \code{NULL}. Random seed for reproducibility. \strong{(both)}
#'   }
#'
#' @return A list with:
#' \itemize{
#'   \item f.predict.hist : n.star \eqn{\times} n.mc matrix. Posterior predictive draws of the latent GP function \eqn{f}.
#'   \item y.predict.hist : n.star \eqn{\times} n.mc matrix. Posterior predictive draws of the response \eqn{y}.
#' }
#'
#' @export
impute.y <- function (Gibbs.data, Gibbs.hist, GP.param, Gibbs.other) {
  # Predict y by kriging (given U.star).
  # single.ustar: whether the goal is to predict at a fixed u.star.
  
  d <- Gibbs.data$d
  y <- Gibbs.data$y
  X <- Gibbs.data$X
  X.star <- Gibbs.data$X.star
  
  sigma.error2.hist <- Gibbs.hist$sigma.error2.hist
  U.hist <-Gibbs.hist$U.hist
  U.star.hist <- Gibbs.hist$U.star.hist
  
  theta <- GP.param$theta
  rho <- GP.param$rho
  
  single.ustar <- if(!is.null(Gibbs.other$single.ustar)) Gibbs.other$single.ustar else FALSE
  separate <- if(!is.null(Gibbs.other$separate)) Gibbs.other$separate else FALSE
  seed <- if(!is.null(Gibbs.other$seed)) Gibbs.other$seed else 222
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Constants ----
  n.mc <- length(sigma.error2.hist)  # number of repeated samples to draw
  n.star <- nrow(X.star)  # number of new samples
  
  # Empty matrices ----
  f.predict.hist <- matrix(0, nrow = n.star, ncol = n.mc)
  y.predict.hist <- matrix(0, nrow = n.star, ncol = n.mc)
  #browser()
  # Sampler ----
  for (idx in 1:n.mc) {
    # Extract ----
    sigma.error2 <- sigma.error2.hist[idx]
    U <- U.hist[, ((idx-1)*d+1):(idx*d), drop = FALSE]
    if (single.ustar) {
      U.star <- U.star.hist
    } else {
      #print(((idx-1)*d+1):(idx*d))
      #browser()
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

