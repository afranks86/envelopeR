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


for(q in 1:qmax) {
    count  <- 1
    for(beta_sd in beta_sd_vec) {
        for(rep in 1:nreps) {      

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

                    envfit1 <- fit_envelope(residual, X, distn="covreg", s=s,
                                            Vinit="OLS",
                                            fmean="YV ~ 1",
                                            fcov="YV ~ X")
                    Ya <- residual
                    
                } else{

                    envfit2 <- fit_envelope(Y, X, distn="covreg", s=s,
                                            fmean = "YV ~ X",
                                            fcov = "YV ~ X",
                                            Vinit="OLS")
                    Ya <- Y
                    
                }

                YV <- Y %*% Vhat
                fit <- covreg.mcmc(YV ~ 1, YV ~ X, R=s)
                cov_psamp <- cov.psamp(fit)
                
                indices  <- sapply(1:nrow(X), function(i) which(apply(unique(X), 1, function(x) all.equal(x, X[i, ]) == "TRUE")))

                cov_psamp <- cov.psamp(envfit2$covariance_list$covreg_res)

                m_psamp  <- covreg::m.psamp(envfit2$covariance_list$covreg_res)

                apply(fit$B1.psamp, 1:2, mean)
                cov_psamp <- cov.psamp(fit)
                beta %*% t(V) %*% Vhat
                
                nsamples <- 100
                SigList <- list()
                for(i in 1:nrow(X)) {
                    SigSamples  <- lapply(1:nsamples, function(s) {
                        Sig  <- cov_psamp[indices[i], , , s]
                    })
                    SigList[[i]]  <- Reduce(`+`, SigSamples)/nsamples
                }

                SigSamples  <- lapply(1:nsamples, function(s) {
                    Sig  <- sum( (cov_psamp[2, , , s] - cov_list[[2]])^2)
                })

                
                envfit <- envfit1
                Vhat <- envfit$V
                Vhat_perp <- rstiefel::NullC(Vhat)

                # Procrustes transformation
                Rsvd <- svd(t(V) %*% Vhat)
                R <- Rsvd$u %*% t(Rsvd$v)
                Vhat <- Vhat %*% t(R)
                
                YVp <- Ya %*% Vhat_perp
                sigma2_hat <- sum(YVp^2)/(n*(p-s))
                
                psi_hat_inv_list <- envfit$covariance_list$SigInvList
                tV <- t(V) %*% Vhat
                proj_cov_hat <- lapply(psi_hat_inv_list, function(psi_inv) {
                    psi <- R %*% solve(psi_inv) %*% t(R)
                    tV %*% psi %*% t(tV) +
                        sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
                })



                
                steins_loss <- sapply(1:n, function(i) {

                    pch <- proj_cov_hat[[i]]
                    Rsvd <- svd(pch %*%  t(cov_list[[i]]))
                    R <- Rsvd$u %*% t(Rsvd$v)
                    
                    
                    steinsLoss(proj_cov_hat[[i]], solve(cov_list[[i]]))
                }) %>% mean

                steins_loss

                subspace_sim  <- tr(envfit$V  %*% t(envfit$V)  %*% V  %*% t(V))/s
                
                subspace_sim_array[a, q, rep, count] <- subspace_sim
                steins_loss_array[a, q, rep, count] <- steins_loss
            }

            count <- count + 1

        }

        save(subspace_sim_array, steins_loss_array, file="sim_results.Rdata")
    }
}
save(subspace_sim_array, steins_loss_array, file="sim_results.Rdata")
