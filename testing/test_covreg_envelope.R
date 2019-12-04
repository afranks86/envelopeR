library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(mvtnorm)
source("../R/fit_envelope.R")
source("../R/envelope_functions.R")



p <- 100
s <- 10
q <- 2
n <- 200

X <- matrix(rnorm(n*q), nrow=n, ncol=q)
    
beta <- matrix(rnorm(q*s, sd=2), nrow=q)

gamma <- matrix(rnorm(q*s, sd=5), nrow=q)


## Covariance: Sigma = gamma %*% XX^T %*%  gamma^T + diag
## Mean = X %*% beta
create_cov <- function(X, gamma) {
    gamma_X <- X %*% gamma
    sig_x <- t(gamma_X) %*%  gamma_X
}

Z <- sapply(1:n, function(i) {
    rmvnorm(1, X[i, ] %*% beta, sigma=create_cov(X[i, ], gamma) + diag(s))
}) %>% t

V  <- rustiefel(p, s)
Y  <- Z %*% t(V)  + rnorm(p * n, mean=0, sd=2)

envfit <- fit_envelope(Y, X, distn="covreg", s=10, Vinit=rustiefel(p, s))

