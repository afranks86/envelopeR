library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(mvtnorm)
library(covreg)
library(envelopeR)


## Number of features
p <- 100

## Rank of the matrix
s <- 4

## Number of predictors
q <- 1
n <- 100

## sd of regression coefficients

create_cov_eigen <- function(X, scaler=1) {

  indices_default  <- 1:3
  Xnew  <- ifelse(indices_default > ncol(X), X[ncol(X)], X[indices_default])

  theta <- pi/2*Xnew[1]
  Lambda <- diag(c(20, 2, 0.5*Xnew[2]+0.5, 0.25*Xnew[3] + 0.25))

  U1 <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
  U2  <- diag(2)
  U  <- Matrix::bdiag(U1, U2)
  sig_X <- U %*% Lambda  %*% t(U)
  as.matrix(sig_X)

}

gamma_sd  <- 1

## error sd
error_sd  <- 0.5


sig_x_rank <- s

beta_sd_vec <- c(0, 1, 3, 5, 8, 10, 20)
nreps <- 1

score_mat  <- matrix(nrow=nreps, ncol=length(beta_sd_vec))


for(i in 1:nreps) {

  X <- matrix(runif(n*q, 0, 1), nrow=n, ncol=q)

  count  <- 1

  for(beta_sd in beta_sd_vec) {

    beta <- matrix(rnorm(q*s, sd=beta_sd), nrow=q)

    ## magnitude of the covariance matrix

    cov_list  <- lapply(1:n, function(i) create_cov_eigen(X[i, , drop=FALSE]))


    Z <- sapply(1:n, function(i) {
      rmvnorm(1, X[i, ] %*% beta, sigma=cov_list[[i]])
    }) %>% t

    V  <- rustiefel(p, s)
    Vnull  <- rstiefel::NullC(V)

    Y  <- Z %*% t(V)  +
      matrix(rnorm(n * (p-s), sd=error_sd), nrow=n, ncol=p-s) %*% t(Vnull)

    envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                           Vinit="OLS")

    subspace_sim  <- tr(envfit$V  %*% t(envfit$V)  %*% V  %*% t(V))/s

    print(subspace_sim)
    score_mat[i, count]  <- subspace_sim
    count  <- count + 1
  }
}
score_mat
