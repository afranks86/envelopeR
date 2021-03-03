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

nbench_times <- 1
Svec <- c(2, 10, 25, 50)
Pvec <- seq(from=100, to=1000, by=100)

## Number of predictors
n <- 100

gamma_sd  <- 1

## error sd
error_sd  <- 1

beta_sd <- 1
nreps <- 10

## Rank of the matrix
q_cur <- 1

create_cov_eigen <- function(X, gammaList, s, scaler=1) {
    
    sig_X <- matrix(0, ncol=s, nrow=s)
    for(i in 1:s) {
        gammaX <- gammaList[[i]] %*% t(X)
        sig_X <- sig_X + tcrossprod(gammaX) 
    }

    sig_X
}

benchmark_array <- array(dim=c(length(Pvec), length(Svec), nbench_times))

pindex <- 1
for(p in Pvec) {
    sindex <- 1
    for(s in Svec) {
        for(rep in 1:nreps) {      

            print("------------------------------")
            print(sprintf("P = %i, S = %i", p, s))
            print("------------------------------")
            
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
                
            benchmarked <- microbenchmark(
                if(s == p) {

                    fit <- covreg.mcmc(Y ~ X - 1, Y ~ X, R=1, verb=FALSE)

                    mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
                    
                    steins_loss <- sapply(1:n, function(i) {
                        a <- t(V) %*% mean_cov_hat[i, , ] %*% V
                        b <- cov_list[[i]]
                        steinsLoss(a, solve(b))
                    }) %>% mean

                    mean_cov_hat <- get_mean_cov_hat(fit)
                    se_loss <- sapply(1:n, function(i) {
                        a <- t(V) %*% mean_cov_hat[i, , ] %*% V
                        b <- cov_list[[i]]
                        sum((a-b)^2)
                    }) %>% mean

                    post_mean <- apply(m.psamp(fit), 1:2, mean)

                    Vhat <- svd(post_mean)$v[, 1:s]
                    
                } else{

                    envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                                           fmean = "YV ~ X - 1",
                                           fcov = "YV ~ X",
                                           Vinit="OLS", tol_v = 1e-4)

                    Vhat <- envfit$V
                    YV <- Y %*% Vhat
                    Vhat_perp <- rstiefel::NullC(Vhat)
                    fit <- covreg.mcmc(YV ~ X - 1, YV ~ X, R=1, verb=FALSE)

                    cov_psamp <- cov.psamp(fit)

                    YVp <- Y %*% Vhat_perp
                    sigma2_hat <- sum(YVp^2)/(n*(p-s))

                },
                times=1)

            ## Benchmarking
            benchmark_array[pindex, sindex, rep] <- benchmarked$time/1e9
        }
        sindex <- sindex + 1
    }
    pindex <- pindex + 1
}

save(benchmark_array, file="benchmark_times.RData")

pdf("paper/Figs/run_times.pdf", height=4, width=8)
dimnames(benchmark_array) = list(S = Svec, P = Pvec, Iter = 1:nbench_times)
as_tibble(reshape2::melt(benchmark_array)) %>% ggplot(aes(x=as.factor(P), y=value)) +
    geom_jitter(aes(col=factor(S))) + scale_fill_discrete_qualitative(palette="Set 2") +
    ylab("Seconds") + ylim(c(0, 300)) +
    xlab("P") + scale_colour_discrete(name  = "S")

dev.off()

dim(benchmark_array)

