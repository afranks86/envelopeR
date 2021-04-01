library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(mvtnorm)

s <- 50
q <- 2
n <- 200

beta_sd  <- 2
gamma_sd  <- 2
error_sd  <- 2

X <- matrix(rnorm(n*q), nrow=n, ncol=q)
beta <- matrix(rnorm(q*s, sd=beta_sd), nrow=q)
gamma <- matrix(rnorm(q*s, sd=gamma_sd), nrow=q)

## Covariance: Sigma = gamma %*% XX^T %*%  gamma^T + diag
## Mean = X %*% beta
create_cov <- function(X, gamma) {
    gamma_X <- X %*% gamma
    sig_x <- t(gamma_X) %*%  gamma_X
}

Y <- sapply(1:n, function(i)
    rmvnorm(1, X[i, ] %*% beta, sigma=create_cov(X[i, ], gamma) + diag(s))) %>% t
svd(Y)$d

cov_dat <- list(X = X,
                 Y = Y,
                 s = s,
                 q = q,
                 n = n)

m1 <- stan_model(file = '../src/stan_files/cov_regression.stan')
m2 <- stan_model(file = '../src/stan_files/cov_regression_aug.stan')

## Fit 1 
nchains <- 4
fit <- stan(file = '../src/stan_files/cov_regression.stan',
            data = cov_dat, chains=nchains,
            init=lapply(1:nchains, function(l)
                list(L_Omega=diag(s),
                     alpha=t(beta),
                     gamma=t(gamma))))

fit <- rstan::vb(m1, data = cov_dat,
                 init=list(L_Omega=diag(s)))

## Bxx^B

colnames(Y)  <- paste0("Y", 1:ncol(Y))
colnames(X)  <- paste0("X", 1:ncol(X))

res  <- covreg::covreg.mcmc(Y ~ X, Y ~ X)

dim(res$B2.psamp)
gsamps  <- res$B2.psamp[, 2:3, 1, ]
plot(t(apply(gsamps, c(1,2), mean)), gamma)

bsamps  <- res$B1.psamp[, 2:3, ]
plot(t(apply(bsamps, c(1,2), mean)), beta)

traceplot(fit, "gamma")

summary(summary(fit)$summary[, "n_eff"])

rstan::get_elapsed_time(fit)

samples <- rstan::extract(fit)

mean((samples$gamma)[, 6, 2])
gamma[2, 6]

apply(samples$gamma, c(2,3), mean)
gamma[1, 1]

plot(t(apply(samples$gamma[1:1000, ,], c(2,3), mean)), gamma)
abline(a=0, b=1)
cor.test(t(apply(samples$gamma, c(2,3), mean)), gamma)



## Fit 1 

nchains <- 4
fit_aug <- stan(file = '../src/stan_files/cov_regression_aug.stan',
            data = cov_dat, chains=nchains,
            init=rep(list(list(L_Omega=diag(s))), nchains),
            iter = 100, warmup = 50)

rstan::get_elapsed_time(fit_aug)

samples <- rstan::extract(fit_aug)

plot(t(apply(samples$gamma, c(2,3), mean)), gamma)

plot(t(apply(samples$beta, c(2,3), mean)), beta)


library(tidyverse)
corrplot::corrplot(cov2cor(apply(samples$A, c(2:3), mean)))
    

stan(file = '../src/stan_files/cov_regression_aug.stan',
            data = cov_dat, chains=nchains,
            init=rep(list(list(L_Omega=diag(s))), nchains),
     iter = 100, warmup = 50)


## Initialize at truth
fit_vb <- rstan::vb(m1, data = cov_dat,
                    init=list(L_Omega=diag(s),
                              alpha=t(beta), gamma=t(gamma)))

## Don't initialize at truth
fit_vb <- rstan::vb(m1, data = cov_dat, iter=20000)

estep  <- covariance_regression_estep(Y, X, m1, method="vb")

mstep  <- covariance_regression_mstep(rustiefel(p, s)

samples <- estep$samples

samples <- rstan::extract(fit_vb)

plot(t(apply(samples$gamma, c(2,3), mean)), gamma)

plot(samples$gamma[1000, , ], t(gamma))
cor.test(t(apply(samples$gamma, c(2,3), mean)), gamma)

plot(t(apply(samples$beta, c(2, 3), mean)), beta)
cor.test(t(apply(samples$beta, c(2, 3), mean)), beta)

fit_vb <- rstan::vb(m1, data = cov_dat,
                    init=list(L_Omega=diag(s)), adapt_engaged=FALSE)


############## GP VB FIT ################

gp_sm <- stan_model(file = '../src/stan_files/gp_cov_regression.stan')
gp_cov_dat <- cov_dat
gp_cov_dat$r <- 2
gp_cov_dat$N <- cov_dat$n
gp_fit <- rstan::vb(gp_sm, data=gp_cov_dat)


gp_samples <- rstan::extract(gp_fit)
dim(gp_samples$beta)

bhat <- t(sapply(1:n, function(i) colMeans(gp_samples$gamma[, i] * gp_samples$beta[, i, ])))
bhat <- apply(gp_samples$beta, c(2, 3), mean)
xg <- X %*% gamma
plot(sapply(1:200, function(i) bhat[i, ] %*% t(bhat[i, ])), 
     sapply(1:200, function(i) xg[i, ] %*% t(xg[i, ])), pch=19, cex=0.2)

cor.test(sapply(1:200, function(i) bhat[i, ] %*% t(bhat[i, ])), 
         sapply(1:200, function(i) xg[i, ] %*% t(xg[i, ])))

gp_samples$rho


fit1 <- rstan::sampling(m1, data=cov_dat, chains=1,
                init=list(list(L_Omega=diag(s))))

fit2 <- rstan::sampling(m2, data=cov_dat, chains=1,
                init=list(list(L_Omega=diag(s))))

fit1_vb <- rstan::vb(m1, data=cov_dat)

samples <- rstan::extract(fit_vb)

plot(t(apply(samples$gamma, c(2,3), mean)), gamma)
cor.test(t(apply(samples$gamma, c(2,3), mean)), gamma)



fit2_vb <- rstan::vb(m2, data=cov_dat, chains=1,
                init=list(list(L_Omega=diag(s))))
