library(tidyverse)
library(magrittr)
library(Matrix)
library(rstiefel)

## This is a general envelope function

## N = sample size
## P = number of features
## S = subspace of mean variation
## R = orthogonal supace of shared covariation
## S + R <= P.  P - S - R dimensions have isotropic covariance sigma^2
## D scales the variances

fit_envelope <- function(Y, X, D = diag(ncol(Y)),
                         s=2, r=0,
                         Vinit = "OLS",
                         Lambda0 = t(X) %*% X,
                         prior_counts=0,
                         Beta0 = matrix(0, nrow=ncol(X), ncol=ncol(Y)),
                         v1=0, v0 = 0,
                         U1=matrix(0, nrow=s, ncol=s),
                         U0=matrix(0, nrow=r, ncol=r),
                         alpha=0, k = 0, center=TRUE, maxIters=1000, use_py=TRUE){

    
  Y <- Y %*% D
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)

  intercept <- rep(0, ncol(Y))
  if(center) {
      intercept = colMeans(Y)
      Y <- scale(Y, scale=FALSE)
  } else {
    intercept <- rep(0, ncol(Y))
  }
  
  Lambda0 <- Lambda0 * prior_counts / n

  if(is.character(Vinit)) {
    if(Vinit == "OLS") {

      beta_hat <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Y
      Vinit <- svd(beta_hat)$v[, 1:min((s+r), q), drop=FALSE]
      if(q < (s+r)) {
        Vinit <- cbind(Vinit, NullC(Vinit)[, 1:(s+r-q), drop=FALSE])
      }
    } else if(Vinit == "COV") {
        
        Vinit <- svd(Y)$v[, 1:min((s+r), q), drop=FALSE]
        if(q < (s+r)) {
            Vinit <- cbind(Vinit, NullC(Vinit)[, 1:(s+r-q), drop=FALSE])
        } 
    } else {
        warn("Randomly initializing")
      Vinit <- rustiefel(ncol(Y), s+r)
    }
  } else if(is.null(Vinit)) {
      Vinit <- rustiefel(ncol(Y), s+r)
  } else if( !is.matrix(Vinit) ) {
      stop("Vinit must be a semi-orthogonal matrix or string")
  } 
  
  if (s <= 0){
    stop("Require s > 0")
  }
  if (r < 0){
    stop("Require r >= 0")
  }
  if(s + r > p) {
    stop("Require s + r <= p")
  }

  if(norm(t(Vinit) %*% Vinit - diag(s+r)) > 1e-6) {
    stop("Vinit not semi-orthogonal")
  }
  
  evals <- eigen(t(X) %*% X)$values
  if (evals[length(evals)] < 1e-6) {
    warning("Perfect Colinearity in X")
  }

  S <- t(Y) %*% Y 
    
  beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
  
  A <- t(Y - X %*% beta_hat) %*% (Y - X %*% beta_hat) + t(Beta0 - beta_hat) %*% Lambda0 %*% (Beta0 - beta_hat)

  ## compute negative log-likelihood
  F <- function(V) {

    G <- V[, 1:s, drop=FALSE]
    if(r > 0) {
      G0 <- V[, (s+1):(s+r)]
      G0part <- (n + r + v0 - 1)/2 * determinant(t(G0) %*% S %*% G0 + U0, logarithm = TRUE)$modulus

      if(r + s < p){
        sig2part <- (n*(p-s-r)/2 + alpha) * log(tr(S)/2 - tr(t(V) %*% S %*% V)/2 + k)
      } else {
        sig2part <- 0
      }
      
    } else {
      G0part <- 0
      sig2part <- (n*(p-s-r)/2 + alpha) * log(tr(S)/2 - tr(t(V) %*% S %*% V)/2 + k)
    }

    ## Minimize the negative log likelihood
    (n + v1 + s +1 - q)/2 * determinant(t(G) %*% A %*% G + U1, logarithm = TRUE)$modulus + G0part + sig2part
                                                                           
  }

  ## compute gradient of negative log-likelihood
  dF <- function(V) {

    G <- V[, 1:s, drop=FALSE]
    Gpart <- - (n + s + v1 + 1 - q) * A %*% G %*% solve(t(G) %*% A %*% G + U1)

    if(r > 0) {
      G0 <- V[, (s+1):(s+r)]
      G0part <- - (n + r + v0 + 1) * S %*% G0 %*% solve(t(G0) %*% S %*% G0 + U0)
    } else {
      G0part <- matrix(0, nrow=nrow(V), ncol=0)
    }

    if(r + s < p){
        sig2part <- - (n*(p-s-r)/2 + alpha) / (tr(S)/2 - tr(t(V) %*% S %*% V)/2 + k) * S %*% V
    } else {
        sig2part <- matrix(0, nrow=p, ncol=s+r)
    }

    

    dV <- cbind(Gpart, matrix(0, nrow=p, ncol=r)) +
      cbind(matrix(0, nrow=p, ncol=s), G0part) +
      sig2part

    ## negative ll grad
    -dV
      
  }

    print("Fitting Stiefel manifold")
    if (use_py) {

        Vfit <- optStiefel_py(Vinit,
                              as.integer(n),
                              as.integer(p),
                              as.integer(s),
                              as.integer(r),
                              as.integer(q),
                              U0, U1, alpha, v0, v1, S, k, A,
                              maxiters = maxIters)
    }
    else {

        Vfit <- optStiefel(
            function(V) F(V),
            function(V) dF(V),
            method = "bb",
            Vinit = Vinit,
            verbose = TRUE,
            maxIters = maxIters,
            maxLineSearchIters = 20
        )
    }

    Yproj <- Y %*% Vfit[, 1:s]
    eta_hat_env <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(Vfit[, 1:s])
    
    list(V=Vfit, intercept=intercept, beta_ols=beta_hat, beta_env=beta_env, eta_hat=eta_hat_env)
    
}


