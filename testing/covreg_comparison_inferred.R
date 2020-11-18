options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(covreg)
library(envelopeR)
library(lubridate)

## Number of features
p <- 25

## Number of predictors
n <- 100

gamma_sd  <- 1

## error sd
error_sd  <- 1

##beta_sd_vec <- 0:5

beta_sd_vec <- c(3)
nreps <- 100

## Rank of the matrix
qmax <- 4

## dimension of v
s <- 8

sfit <- c(2:s)

create_cov_eigen <- function(X, gammaList, s, scaler=1) {
    
    sig_X <- matrix(0, ncol=s, nrow=s)
    for(i in 1:s) {
        gammaX <- gammaList[[i]] %*% t(X)
        sig_X <- sig_X + tcrossprod(gammaX) 
    }

    sig_X
}

file_name <- sprintf("sim_s_comparison_inferred_%s.Rdata", lubridate::today())
arr_dimnames <- list("s" = 1:3, "type" = c("s_hat", "s", "p"), "rep" = 1:nreps)
steins_loss_array  <- array(dim=c(length(sfit), 3, nreps), dimnames=arr_dimnames)
squared_error_loss_array  <- array(dim=c(length(sfit), 3, nreps), dimnames=arr_dimnames)

for(rep in 1:nreps) {      

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

    s_count <- 1
    
    for(s_cur in sfit) {

        ## s = s_cur
        envfit <- fit_envelope(Y, X, distn="covreg", s=s_cur,
                               fmean = "YV ~ X - 1",
                               fcov = "YV ~ X",
                               Vinit="OLS")

        Vstar <- envfit$V
        YV <- Y %*% Vstar
        Vstar_perp <- rstiefel::NullC(Vstar)
        fit <- covreg.mcmc(YV ~ X - 1, YV ~ X, R=s, verb=FALSE)

        cov_psamp <- cov.psamp(fit)

        YVp <- Y %*% Vstar_perp
        sigma2_hat <- sum(YVp^2)/(n*(p-s))

        mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
        steins_loss <- sapply(1:n, function(i) {
            b <- t(Vstar) %*%  (V %*% cov_list[[i]] %*% t(V) + error_sd^2 * Vnull %*% t(Vnull))%*% Vstar
            steinsLoss(mean_cov_hat[i, ,], solve(b))
        }) %>% mean

        mean_cov_hat <- get_mean_cov_hat(fit)
        se_loss <- sapply(1:n, function(i) {
            b <- t(Vstar) %*%  (V %*% cov_list[[i]] %*% t(V) + error_sd^2 * Vnull %*% t(Vnull))%*% Vstar
            sum((mean_cov_hat[i, , ] - b)^2)
        }) %>% mean


        steins_loss_array[s_count, 1, rep] <- steins_loss
        squared_error_loss_array[s_count, 1, rep] <- se_loss
        
        ## s = s
        envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                               fmean = "YV ~ X - 1",
                               fcov = "YV ~ X",
                               Vinit="OLS")


        Vhat <- envfit$V
        YV <- Y %*% Vhat
        Vhat_perp <- rstiefel::NullC(Vhat)
        fit <- covreg.mcmc(YV ~ X - 1, YV ~ X, R=s, verb=FALSE)

        cov_psamp <- cov.psamp(fit)

        YVp <- Y %*% Vhat_perp
        sigma2_hat <- sum(YVp^2)/(n*(p-s))

        mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
        steins_loss <- sapply(1:n, function(i) {
            a <- t(Vstar) %*% (Vhat %*% mean_cov_hat[i, ,] %*% t(Vhat) + sigma2_hat * Vhat_perp %*% t(Vhat_perp)) %*% Vstar  
            b <- t(Vstar) %*%  (V %*% cov_list[[i]] %*% t(V) + error_sd^2 * Vnull %*% t(Vnull))%*% Vstar
            steinsLoss(a, solve(b))
        }) %>% mean

        mean_cov_hat <- get_mean_cov_hat(fit)
        se_loss <- sapply(1:n, function(i) {
            a <- t(Vstar) %*% (Vhat %*% mean_cov_hat[i, ,] %*% t(Vhat) + sigma2_hat * Vhat_perp %*% t(Vhat_perp)) %*% Vstar  
            b <- t(Vstar) %*%  (V %*% cov_list[[i]] %*% t(V) + error_sd^2 * Vnull %*% t(Vnull))%*% Vstar
            sum((a-b)^2)
        }) %>% mean

        steins_loss_array[s_count, 2, rep] <- steins_loss
        squared_error_loss_array[s_count, 2, rep] <- se_loss
        
        ## s = p
        fit <- covreg.mcmc(Y ~ X - 1, Y ~ X, R=p, verb=FALSE)

        post_mean <- apply(m.psamp(fit), 1:2, mean)

        Vhat <- svd(post_mean)$v[, 1:s]

        mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
        steins_loss <- sapply(1:n, function(i) {
            a <- t(Vstar) %*% mean_cov_hat[i, ,] %*% Vstar
            b <- t(Vstar) %*%  (V %*% cov_list[[i]] %*%  t(V) + error_sd^2 * Vnull %*% t(Vnull)) %*% Vstar
            steinsLoss(a, solve(b))
        }) %>% mean

        mean_cov_hat <- get_mean_cov_hat(fit)
        se_loss <- sapply(1:n, function(i) {
            a <- t(Vstar) %*% mean_cov_hat[i, ,] %*% Vstar
            b <- t(Vstar) %*%  (V %*% cov_list[[i]] %*%  t(V) + error_sd^2 * Vnull %*% t(Vnull)) %*% Vstar
            sum((a-b)^2)
        }) %>% mean

        steins_loss_array[s_count, 3, rep] <- steins_loss
        squared_error_loss_array[s_count, 3, rep] <- se_loss
        
        ## 
        print(sprintf("--------- %i, %i ----------------", s_cur, rep))

        s_count <- s_count +1
        
    }

    save(steins_loss_array, squared_error_loss_array, file=file_name)
}


