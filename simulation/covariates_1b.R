###################################
####### Figure 2b in Franks 2020
##################################

options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(covreg)
library(envelopeR)


create_cov_eigen <- function(X, scaler=1) {

  indices_default  <- 1:3
  Xnew  <- ifelse(indices_default > ncol(X), X[ncol(X)], X[indices_default])

  theta <- pi/2*Xnew[1]
  Lambda <- diag(c(20, 2, 0.5*Xnew[2]+0.5, 0.25*Xnew[3] + 0.25))*scaler

  U1 <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
  U2  <- diag(2)
  U  <- as.matrix(Matrix::bdiag(U1, U2))
  sig_X <- U %*% Lambda  %*% t(U)
  as.matrix(sig_X)

}

## Number of features
p <- 100

## Number of predictors
n <- 100

gamma_sd  <- 1

## error sd
error_sd  <- 0.5

beta_sd_vec <- 0:10
nreps <- 100

## Rank of the matrix
qmax <- 4

## dimension of v
s <- 4

score_array  <- array(dim=c(s, nreps, length(beta_sd_vec)))

for(q in 1:4) {
    for(i in 1:nreps) {

        X <- matrix(runif(n*q, 0, 1), nrow=n, ncol=q)
        cov_list  <- lapply(1:n, function(i) create_cov_eigen(X[i, , drop=FALSE], scaler=1))
        V  <- rustiefel(p, s)
        Vnull  <- rstiefel::NullC(V)

        beta_mat <- matrix(runif(q*s, 1, 2), nrow=q)

        count  <- 1

        for(beta_sd in beta_sd_vec) {

            beta <- beta_sd * beta_mat

            Z <- sapply(1:n, function(i) {
                rmvnorm(1, X[i, ] %*% beta, sigma=cov_list[[i]])
            }) %>% t

            Y  <- Z %*% t(V)  +
                matrix(rnorm(n * (p-s), sd=error_sd), nrow=n, ncol=p-s) %*% t(Vnull)

            envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                                   Vinit="OLS")

            subspace_sim  <- tr(envfit$V  %*% t(envfit$V)  %*% V  %*% t(V))/s

            score_array[q, i, count]  <- subspace_sim
            count  <- count + 1
        }
    }

    save(score_array, file="scores.Rdata")
}

save(score_array, file="scores.Rdata")
