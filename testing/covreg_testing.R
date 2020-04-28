options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(covreg)
library(envelopeR)


create_cov_eigen <- function(X, gammaList, s, scaler=1) {
    
    sig_X <- matrix(0, ncol=s, nrow=s)
    for(i in 1:s) {
        gammaX <- gammaList[[i]] %*% t(X)
        sig_X <- sig_X + tcrossprod(gammaX)
    }

    sig_X
}


get_mean_cov_hat <- function(fit, inv=FALSE) 
{
    A.psamp = fit$A.psamp
    B.psamp = fit$B2.psamp
    nsave = dim(A.psamp)[3]
    p = dim(A.psamp)[1]
    R = dim(B.psamp)[3]
    X = unique(fit$matrix.cov)
    n = dim(X)[1]
    s.psamp = array(0, dim = c(n, p, p))
    for (i in 1:n) {
        for (iter in 1:nsave) {
            ss = A.psamp[, , iter]
            for (r in 1:R) {
                ss = ss + B.psamp[, , r, iter] %*% X[i, ] %*% 
                  t(X[i, ]) %*% t(B.psamp[, , r, iter])
            }
            if(inv)
                s.psamp[i, , ] = s.psamp[i, , ] + solve(ss)
            else
                s.psamp[i, , ] = s.psamp[i, , ] + ss
        }
        if(inv)
            s.psamp[i, , ] <- solve(s.psamp[i, , ]/nsave)
        else
            s.psamp[i, , ] <- s.psamp[i, , ]/nsave
    }
    s.psamp

}


## Number of features
p <- 100

## Number of predictors
n <- 100

gamma_sd  <- 1

## error sd
error_sd  <- .1

beta_sd_vec <- seq(2, 10, by=2)

nreps <- 100

## Rank of the matrix
qmax <- 1

## dimension of v
s <- 4

subspace_sim_array  <- array(dim=c(2, s, nreps, length(beta_sd_vec)))
steins_loss_array  <- array(dim=c(2, s, nreps, length(beta_sd_vec)))
squared_error_loss_array  <- array(dim=c(2, s, nreps, length(beta_sd_vec)))
for(rep in 1:nreps) {      
    for(q in 1:qmax) {
        count  <- 1
        for(beta_sd in beta_sd_vec) {

            X <- matrix(rnorm(n*q, 0, 1), nrow=n, ncol=q)

            gammaList <- lapply(1:s, function(i) matrix(rnorm(s*q, 0, sd=gamma_sd), nrow=s, ncol=q))
            
            cov_list  <- lapply(1:n, function(i) create_cov_eigen(X[i, , drop=FALSE], gammaList, s, scaler=1) + error_sd^2*diag(s))
            V  <- rustiefel(p, s)
            Vnull  <- rstiefel::NullC(V)

            beta_mat <- matrix(rnorm(q*s, 0, 1), nrow=q)
            beta <- beta_sd * beta_mat

            Z <- sapply(1:n, function(i) {
                rmvnorm(1, X[i, ] %*% beta, sigma=cov_list[[i]])
            }) %>% t
            
            Y  <- Z %*% t(V)  +
                matrix(rnorm(n * (p-s), sd=10*error_sd), nrow=n, ncol=p-s) %*% t(Vnull)
            
            for(a in 1:2) {

                if(a == 1) {
                    residual <- matrix(0, nrow=n, ncol=p)
                    for(i in 1:p)
                        residual[, i] <- lm(Y[, i] ~ X)$residuals

                    envfit <- fit_envelope(residual, X, distn="covreg", s=s,
                                            Vinit="OLS",
                                            fmean="YV ~ 1",
                                            fcov="YV ~ X")
                    Ya <- residual
                    Vhat <- envfit$V
                    YV <- Ya %*% Vhat
                    Vhat_perp <- rstiefel::NullC(Vhat)
                    fit <- covreg.mcmc(YV ~ 1, YV ~ X, R=s, verb=FALSE)
                    
                } else{

                    envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                                            fmean = "YV ~ X - 1",
                                            fcov = "YV ~ X",
                                           Vinit="OLS")
                    Ya <- Y

                    Vhat <- envfit$V
                    YV <- Ya %*% Vhat
                    Vhat_perp <- rstiefel::NullC(Vhat)
                    fit <- covreg.mcmc(YV ~ X - 1, YV ~ X, R=s, verb=FALSE)
                }
                
                m_psamp  <- covreg::m.psamp(fit)

                YVp <- Ya %*% Vhat_perp
                sigma2_hat <- sum(YVp^2)/(n*(p-s))

                mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
                steins_loss <- sapply(1:n, function(i) {
                    a <- t(V) %*% Vhat %*% mean_cov_hat[i, ,] %*% t(Vhat) %*% V  + sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
                    b <- cov_list[[i]]
                    steinsLoss(a, solve(b))
                }) %>% mean

                mean_cov_hat <- get_mean_cov_hat(fit)
                se_loss <- sapply(1:n, function(i) {
                    a <- t(V) %*% Vhat %*% mean_cov_hat[i, ,] %*% t(Vhat) %*% V  + sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
                    b <- cov_list[[i]]
                    sum((a-b)^2)
                }) %>% mean


                subspace_sim  <- tr(Vhat  %*% t(Vhat)  %*% V  %*% t(V))/s
                
                subspace_sim_array[a, q, rep, count] <- subspace_sim
                steins_loss_array[a, q, rep, count] <- steins_loss
                squared_error_loss_array[a, q, rep, count] <- se_loss
            }


            print(sprintf("-------------------------------------------------- No Mean:, %i, %i: %f, %f, %f ----------------", rep, beta_sd, subspace_sim_array[1, q, rep, count], steins_loss_array[1, q, rep, count], squared_error_loss_array[1, q, rep, count]))
            print(sprintf("--------------------------------------------------- With Mean:, %i, %i: %f, %f, %f ----------------", rep, beta_sd, subspace_sim, steins_loss, se_loss))
            count <- count + 1
        }

        save(subspace_sim_array, steins_loss_array, squared_error_loss_array, file="sim_mean_results.Rdata")
    }
}
save(subspace_sim_array, steins_loss_array, squared_error_loss_array, file="sim_mean_results.Rdata")
