rm(list=ls())

options(mc.cores = parallel::detectCores())

## Need these for Rscript
library(microbenchmark)
library(rstiefel)
library(tidyverse)
library(cowplot)
library(colorspace)
library(mvtnorm)
library(glmnet)
library(covreg)

devtools::build("../../envelopeR")
devtools::install("../../envelopeR")

library(envelopeR)

s <- 4
p <- 200

## 1, 1, 3 was all about the same
## gamma_sd = 0.1, error_sd=1, beta_sd=1 OLS does best, cov second best
## gamma_sd = 0.75, error_sd=1, beta_sd=.25 ...

## Number of predictors
n <- 100

gamma_sd  <- 1

## error sd
error_sd  <- 1

nreps <- 10

## Rank of the matrix
q_cur <- 4

create_cov_eigen <- function(X, gammaList, s, scaler=1) {
    
    sig_X <- matrix(0, ncol=s, nrow=s)
    for(i in 1:s) {
        gammaX <- gammaList[[i]] %*% t(X)
        sig_X <- sig_X + tcrossprod(gammaX) 
    }

    sig_X
}

beta_vec <- c(0.1, 0.5, 0.75, 0.75, 0.5,, 2)
gamma_vec <- c(0.5, 0.5, 0.1, 0.75, 1, 1)

params <- mapply(function(x, y) c(beta=x, gamma=y), beta_vec, gamma_vec, SIMPLIFY=FALSE)

subspace_init_mat <- subspace_sim_mat <- steins_loss_mat <- array(dim=c(length(params), 2, nreps))


for(rep in 1:nreps) {
    count <- 0
    for(par in params){
        count <- count + 1

        beta_sd <- par[["beta"]]
        gamma_sd <- par[["gamma"]]


        X <- matrix(rnorm(n*q_cur, 0, 1), nrow=n, ncol=q_cur)

        gammaList <- lapply(1:s, function(i) matrix(rnorm(s*q_cur, 0, sd=gamma_sd), nrow=s, ncol=q_cur))
        cov_list  <- lapply(1:n, function(i) create_cov_eigen(X[i, , drop=FALSE], gammaList, s, scaler=1) + error_sd^2 * diag(s))

        V  <- rustiefel(p, s)
        Vnull  <- rstiefel::NullC(V)

        beta_mat <- matrix(rnorm(q_cur*s, 0, 1), nrow=q_cur)
        beta <- beta_sd * beta_mat

        Z <- sapply(1:n, function(i) {
            rmvnorm(1, X[i, ] %*% beta, sigma=cov_list[[i]])
        }) %>% t

        Y  <- Z %*% t(V)  +
            matrix(rnorm(n * (p-s), sd=error_sd), nrow=n, ncol=p-s) %*% t(Vnull)

        beta_hat <- solve(t(X) %*% X) %*% (t(X) %*% Y)
        resid <- Y - X %*% beta_hat

        
        for(type in 1:2) {
            print("------------------------------")
            print(sprintf("Type = %i, Rep = %i", type, rep))
            print("------------------------------")

            if(type==1){ 
                Vinit <- rustiefel(p, s)
                subspace_init_mat[count, 1, rep] <- tr(Vinit  %*% t(Vinit)  %*% V  %*% t(V))/s
            } else  {
                Vinit <- svd(beta_hat)$v[, 1:s, drop=FALSE]
                subspace_init_mat[count, 2, rep] <- tr(Vinit  %*% t(Vinit)  %*% V  %*% t(V))/s
            }
            
            envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                                   fmean = "YV ~ X - 1",
                                   fcov = "YV ~ X",
                                   Vinit=Vinit, tol_v = 1e-4, nchunks=1)

            Vhat <- envfit$V
            YV <- Y %*% Vhat
            Vhat_perp <- rstiefel::NullC(Vhat)
            fit <- covreg.mcmc(YV ~ X - 1, YV ~ X, R=1, verb=FALSE)

            cov_psamp <- cov.psamp(fit)

            YVp <- Y %*% Vhat_perp
            sigma2_hat <- sum(YVp^2)/(n*(p-s))

            mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
            steins_loss <- sapply(1:n, function(i) {
                a <- t(V) %*% Vhat %*% mean_cov_hat[i, ,] %*% t(Vhat) %*% V  + sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
                b <- cov_list[[i]]
                steinsLoss(a, solve(b))
            }) %>% mean

            subspace_sim  <- tr(Vhat  %*% t(Vhat)  %*% V  %*% t(V))/s
            subspace_sim
            
            subspace_sim_mat[count, type, rep] <- subspace_sim
            steins_loss_mat[count, type, rep] <- steins_loss
            print(subspace_sim)
        }
    }
}

save(subspace_sim_mat, steins_loss_mat, file="init_sims_full.RData")

subspace_init_mat
subspace_sim_mat
## steins_loss_mat
## steinsLoss


## summary(t(subspace_sim_mat))

## Vinit <- rustiefel(p, s)
## print(tr(Vinit  %*% t(Vinit)  %*% V  %*% t(V))/s)
## envfit <- fit_envelope(Y, X, distn="covreg", s=s,
##                        fmean = "YV ~ X - 1",
##                        fcov = "YV ~ X",
##                        Vinit=rustiefel(p, s), tol_v = 1e-4)
## envfit <- fit_envelope(Y, X, distn="covreg", s=s,
##                        fmean = "YV ~ X - 1",
##                        fcov = "YV ~ X",
##                        Vinit="OLS", tol_v = 1e-4)

## Vhat <- envfit$V
## tr(Vhat  %*% t(Vhat)  %*% V  %*% t(V))/s

## Vinit <- svd(beta_hat)$v[, 1:s, drop=FALSE]

## Vinit <- svd(Y)$v[, 1:s, drop=FALSE]
## tr(Vinit  %*% t(Vinit)  %*% V  %*% t(V))/s
