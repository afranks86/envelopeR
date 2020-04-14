options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(covreg)
library(envelopeR)


## create_cov_eigen <- function(X, scaler=1) {

##     indices_default  <- 1:3
##     Xnew  <- ifelse(indices_default > ncol(X), X[ncol(X)], X[indices_default])

##     theta <- pi/2*Xnew[1]
##     Lambda <- diag(c(20, 2, 0.5*Xnew[2]+0.5, 0.25*Xnew[3] + 0.25))*scaler

##     U1 <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
##     U2  <- diag(2)
##     U  <- as.matrix(Matrix::bdiag(U1, U2))
##     sig_X <- U %*% Lambda  %*% t(U)
##     as.matrix(sig_X)

## }


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

subspace_sim_array  <- array(dim=c(2, s, nreps, length(beta_sd_vec)))
steins_loss_array  <- array(dim=c(2, s, nreps, length(beta_sd_vec)))
squared_error_loss_array  <- array(dim=c(2, s, nreps, length(beta_sd_vec)))
for(rep in 1:nreps) {      
    for(q in 1:qmax) {
        count  <- 1
        for(beta_sd in beta_sd_vec) {

            X <- matrix(runif(n*q, 0, 1), nrow=n, ncol=q)
            cov_list  <- lapply(1:n, function(i) create_cov_eigen(X[i, , drop=FALSE], scaler=1))
            V  <- rustiefel(p, s)
            Vnull  <- rstiefel::NullC(V)

            beta_mat <- matrix(runif(q*s, 1, 2), nrow=q)

            count  <- 1

            beta <- beta_sd * beta_mat

            Z <- sapply(1:n, function(i) {
                rmvnorm(1, X[i, ] %*% beta, sigma=cov_list[[i]])
            }) %>% t
            
            Y  <- Z %*% t(V)  +
                matrix(rnorm(n * (p-s), sd=error_sd), nrow=n, ncol=p-s) %*% t(Vnull)
            
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
                    fit <- covreg.mcmc(YV ~ 1, YV ~ X, R=s)
                    
                } else{

                    envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                                            fmean = "YV ~ X - 1",
                                            fcov = "YV ~ X",
                                            Vinit="OLS")
                    Ya <- Y

                    Vhat <- envfit$V
                    YV <- Ya %*% Vhat
                    Vhat_perp <- rstiefel::NullC(Vhat)
                    fit <- covreg.mcmc(YV ~ X - 1, YV ~ X, R=s)
                }
                
                
                cov_psamp <- cov.psamp(fit)
                m_psamp  <- covreg::m.psamp(fit)

                YVp <- Ya %*% Vhat_perp
                sigma2_hat <- sum(YVp^2)/(n*(p-s))
                
                steins_loss <- sapply(1:n, function(i) {
                    a <- t(V) %*% Vhat %*% apply(cov_psamp, 1:3, mean)[i, ,] %*% t(Vhat) %*% V  + sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
                    b <- cov_list[[i]]
                    steinsLoss(a, solve(b))
                }) %>% mean

                se_loss <- sapply(1:n, function(i) {
                    a <- t(V) %*% Vhat %*% apply(cov_psamp, 1:3, mean)[i, ,] %*% t(Vhat) %*% V  + sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
                    b <- cov_list[[i]]
                    sum((a-b)^2)
                }) %>% mean


                subspace_sim  <- tr(Vhat  %*% t(Vhat)  %*% V  %*% t(V))/s
                
                subspace_sim_array[a, q, rep, count] <- subspace_sim
                steins_loss_array[a, q, rep, count] <- steins_loss
                squared_error_loss_array[a, q, rep, count] <- se_loss
            }

            count <- count + 1

        }

        save(subspace_sim_array, steins_loss_array, squared_error_loss_array, file="sim_results2.Rdata")
    }
}
save(subspace_sim_array, steins_loss_array, squared_error_loss_array, file="sim_results2.Rdata")
