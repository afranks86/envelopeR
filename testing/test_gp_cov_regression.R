library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(mvtnorm)

s <- 3
q <- 2
n <- 100

X <- matrix(rnorm(n*q), nrow=n, ncol=q)
    
eta <- matrix(rnorm(q*s, sd=2), nrow=q)

beta <- matrix(rnorm(q*s, sd=2), nrow=q)


get_cov <- function(X, beta) {
    beta_X <- X %*% beta
    sig_x <- t(beta_X) %*%  beta_X
}

Y <- sapply(1:n, function(i)
    rmvnorm(1, X[i, ] %*% eta, sigma=get_cov(X[i, ], beta) + diag(s))) %>% t
svd(Y)$d

cov_dat <- list(X = X,
                Y = Y,
                s = s,
                q = q,
                N = n,
                r = 1 )

m1 <- stan_model(file = '../src/stan_files/gp_cov_regression.stan')

fit1 <- rstan::sampling(m1, data=cov_dat, chains=1)

plot(fit1)
samples <- rstan::extract(fit1)
dim(samples$beta)



B <- apply(samples$beta, 2:3, mean)


idx <- 1
BX <- t(B[idx, , drop=FALSE]) %*% B[idx, , drop=FALSE]
plot(t(X[1, ] %*% beta) %*% (X[1, ] %*%  beta), t(BX) %*% BX)

abline(a=0, b=1)


## Fit 1 
nchains <- 4
fit <- stan(file = '../src/stan_files/cov_regression.stan',
            data = cov_dat, chains=nchains,
            init=rep(list(list(L_Omega=diag(s))), nchains),
            iter = 100, warmup = 50)

rstan::get_elapsed_time(fit)

samples <- rstan::extract(fit)

plot(t(apply(samples$gamma, c(2,3), mean)), beta)

plot(t(apply(samples$beta, c(2,3), mean)), eta)


library(rstan)
library(shinystan)
library(spNNGP)       # Build neighbor index
source("../../nngp/NNGP_STAN/NNmatrix.R")  # Build matrix including nearest neighbor information
nngp_mod <- stan_model(file = '../src/stan_files/nngp_cov_regression.stan')


#---------------------- Build neighbor index by NNMatrix ----------------------#
M = 6                 # Number of Nearest Neighbors


Xtilde <- cbind(X, X)
NN.matrix <- NNMatrix(coords = Xtilde[, 1:2], n.neighbors = M, n.omp.threads = 2)
str(NN.matrix)
Check_Neighbors(NN.matrix$coords.ord, n.neighbors = M, NN.matrix, ind = 100)

#-------------------------- Set parameters of priors --------------------------#

ss = 3 * sqrt(2)       # scale parameter in the normal prior of sigma 
st = 3 * sqrt(0.1)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       # shape and rate parameters in the Gamma prior of phi 

#------------------------------ NNGP response ---------------------------------#





cov_dat <- list(X = X[NN.matrix$ord, 1:q, drop=FALSE],
                Y = Y[NN.matrix$ord, ],
                s = s,
                q = q,
                N = n,
                r = 1,
                M = M,
                NN_ind = NN.matrix$NN_ind,
                NN_dist = NN.matrix$NN_dist,
                NN_distM = NN.matrix$NN_distM,
                ss = ss, st = st, ap = ap, bp = bp)

myinits <-list(list(sigma = 10, tau = 0.01, phi = 12),
               list(sigma = 15, tau = 0.01, phi = 5),
               list(sigma = 25, tau = 0.01, phi = 9))


nngp_fit <- rstan::sampling(nngp_mod, data=cov_dat, init=myinits,
                            iter = 1000, chains=3)

nngp_samples <- rstan::extract(nngp_fit)

plot(apply(nngp_samples$beta, c(2, 3), mean)[, 20],
     (X[NN.matrix$ord, , drop=FALSE] %*% beta)[, 20])

cor.test(apply(nngp_samples$beta, c(2, 3), mean),
         X[NN.matrix$ord, , drop=FALSE] %*% beta)


bhat <- apply(nngp_samples$beta, c(2, 3), mean)
xb <- X[NN.matrix$ord, , drop=FALSE] %*% beta
plot(sapply(1:200, function(i) bhat[i, ] %*% t(bhat[i, ])), 
     sapply(1:200, function(i) xb[i, ] %*% t(xb[i, ])))

plot(bhat %*% t(bhat), 


mean(nngp_samples$sigmasq)




plot(apply(nngp_samples$beta, c(2, 3), mean)[, 1])
plot((X[NN.matrix$ord, ,drop=FALSE] %*% beta)[, 1])


