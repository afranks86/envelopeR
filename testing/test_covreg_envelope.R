library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(mvtnorm)
library(covreg)

source("R/fit_envelope.R")
source("R/envelope_functions.R")


## Number of features
p <- 50

## Rank of the matrix
s <- 5

## Number of predictors
q <- 3
n <- 1000

## sd of regression coefficients
beta_sd  <- 0.2

## magnitude of the covariance matrix
gamma_sd  <- 2

## error sd
error_sd  <- 1

X <- matrix(rnorm(n*q), nrow=n, ncol=q)

## Test more observatiosn near 0?

X <- matrix(rnorm(n/2*q), nrow=n/2, ncol=q)
X <- rbind(X, matrix(runif(n/2*q, -0.1, 0.1), nrow=n/2, ncol=q))


X <- matrix(runif(n*q, 0, 2), nrow=n, ncol=q)

## Regression coefficients
beta <- matrix(rnorm(q*s, sd=beta_sd), nrow=q)

## Covariance regression coefficients
gamma <- matrix(rnorm(q*s, sd=gamma_sd), nrow=q)


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
Vnull  <- rstiefel::NullC(V)

Y  <- Z %*% t(V)  +
    matrix(rnorm(n * (p-s), sd=error_sd), nrow=n, ncol=p-s) %*% t(Vnull)

## Random initialization
envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                       Vinit=rustiefel(p, s))
tr(envfit$V %*% t(envfit$V) %*% V %*% t(V)) / s
envfit_normal <- fit_envelope(Y, X, distn="normal", s=s,
                              Vinit=rustiefel(p, s))
tr(envfit_normal$V %*% t(envfit_normal$V) %*% V %*% t(V)) / s

## OLS and/or COV initialization (normal)
envfit_normal <- fit_envelope(Y, X, distn="normal", s=s,
                              Vinit="OLS")
tr(envfit_normal$V %*% t(envfit_normal$V) %*% V %*% t(V)) / s
envfit_normal <- fit_envelope(Y, X, distn="normal", s=s,
                              Vinit="COV")
tr(envfit_normal$V %*% t(envfit_normal$V) %*% V %*% t(V)) / s


## OLS and/or COV initialization (cov-reg)
envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                       Vinit="OLS")
tr(envfit$V %*% t(envfit$V) %*% V %*% t(V)) / s
envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                       Vinit="COV")
tr(envfit$V %*% t(envfit$V) %*% V %*% t(V)) / s


envfit_normal <- fit_envelope(Y, X, distn="normal", s=s,
                       Vinit=envfit$V)

envfit$F(envfit$V)
envfit_normal$F(envfit$V)
envfit_normal$F(envfit_normal$V)

sum((envfit$beta_env - beta %*% t(V))^2)
sum((envfit$beta_ols - beta %*% t(V))^2)



estep  <- covariance_regression_estep(Y %*% envfit$V, X,
                                      method="covreg")
estep <- envfit$covariance_list

### TO DO: recover gammas?
dim(estep$samples$gamma)
ghat  <- apply(estep$samples$gamma, c(1, 2), mean)

covx  <- gamma  %*% t(gamma)
covx_hat  <- t(ghat)  %*% ghat

plot(covx, covx_hat)
cor.test(covx, covx_hat)

names(estep$samples)


lapply(1:n, function(j) {
    j_lst  <- lapply(1:1000, function(i) {
        gx  <- estep$samples$gamma[, , i] * X[j, ]
        gx %*% t(gx)
    })
    Sig_x <- Reduce(j_lst, '+')/1000
    Sig_x
})
plot(apply(abs(estep$samples$gamma), c(2, 3), mean), abs(gamma))

