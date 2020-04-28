library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(mvtnorm)
library(covreg)
library(envelopeR)


## Number of features
p <- 100

## Rank of the matrix
s <- 4

## Number of predictors
q <- 1
n <- 100

## sd of regression coefficients


beta_sd  <- 0


## magnitude of the covariance matrix
gamma_sd  <- 1

## error sd
error_sd  <- 0.5

X <- matrix(rnorm(n*q), nrow=n, ncol=q)

## ## Test more observations near 0?
## X <- matrix(rnorm(n/2*q), nrow=n/2, ncol=q)
## X <- rbind(X, matrix(runif(n/2*q, -0.1, 0.1), nrow=n/2, ncol=q))


X <- matrix(runif(n*q, 0, 1), nrow=n, ncol=q)

## Regression coefficients
beta <- matrix(rnorm(q*s, sd=beta_sd), nrow=q)

## Covariance regression coefficients

sig_x_rank <- s

gamma_list <- lapply(1:sig_x_rank, function(i) matrix(rnorm(q*s, sd=gamma_sd), nrow=q))
gamma_intercept_list  <- lapply(1:sig_x_rank, function(i) rnorm(s, sd=2*gamma_sd))

## Covariance: Sigma = gamma %*% XX^T %*%  gamma^T + diag
## Mean = X %*% beta
## gamma_i is for intercept
create_cov <- function(X, gamma_list, gamma_intercept_list=NULL) {

  if(is.null(gamma_intercept))
    sig_X <- Reduce('+', lapply(1:sig_x_rank, function(i) {
      gamma_x <- X %*% gamma_list[[i]]
      sig_x <- t(gamma_X) %*%  gamma_X
      sig_x
    }))
  else
    sig_X <- Reduce('+', lapply(1:sig_x_rank, function(i) {
      gamma_x <- c(1, X)  %*% rbind(gamma_intercept_list[[i]], gamma_list[[i]])
      sig_x <- t(gamma_x) %*%  gamma_x
      sig_x
    }))

  sig_X

}

create_cov_eigen <- function(X, scaler=1) {

  indices_default  <- 1:3
  Xnew  <- ifelse(indices_default > ncol(X), X[ncol(X)], X[indices_default])

  theta <- pi/2*Xnew[1]
  Lambda <- diag(c(200, 20, 5*Xnew[2]+5, 2.5*Xnew[3] + 2.5)) * scaler

  U1 <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
  U2  <- diag(2)
  U  <- Matrix::bdiag(U1, U2)
  sig_X <- U %*% Lambda  %*% t(U)
  as.matrix(sig_X)

}


cov_list  <- lapply(1:n, function(i) create_cov(X[i, , drop=FALSE], gamma_list, gamma_intercept_list))

cov_list  <- lapply(1:n, function(i) create_cov_eigen(X[i, , drop=FALSE], scaler=1))


Z <- sapply(1:n, function(i) {
    rmvnorm(1, X[i, ] %*% beta, sigma=cov_list[[i]])
}) %>% t

V  <- rustiefel(p, s)
Vnull  <- rstiefel::NullC(V)

Y  <- Z %*% t(V)  +
    matrix(rnorm(n * (p-s), sd=error_sd), nrow=n, ncol=p-s) %*% t(Vnull)


envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                       Vinit="COV")

tr(envfit$V  %*% t(envfit$V)  %*% V  %*% t(V))/s


cov_psamp  <- covreg::cov.psamp(envfit$covariance_list$covreg_res)


Vstar <- envfit$V[, 1:2]
O <- t(envfit$V)  %*% Vstar

nsamps  <- dim(cov_psamp)[4]

Osamps_proj <- array(dim = c(2, 2, n, nsamps))
omegaSamps_proj <- array(dim = c(2, n, nsamps))
cov_proj <- array(dim = c(n, 2, 2, nsamps))

to_plot <- 1:n
for(i in 1:nsamps) {
  for(k in 1:length(to_plot)) {

    cov_proj_ik <- t(O) %*% cov_psamp[to_plot[k], , , i]  %*% O
    cov_proj[k, , , i]  <- cov_proj_ik
    eig <- eigen(cov_proj_ik)
    Osamps_proj[, , k, i] <- eig$vectors
    lambda <- eig$values
    omegaSamps_proj[, k, i] <- lambda/(lambda+1)
  }
}


## obs_to_plot <- sample(100, 3)
ix <- sort(X[, 1], index.return=TRUE)$ix
obs_to_plot <- ix[c(1, 50, 100)]
posterior_plot <- posteriorPlot(cov_proj,
                                Osamps_proj, omegaSamps_proj,
                                nsamps=nsamps,
                                obs_to_plot=obs_to_plot,
                                type="logratio",
                                ymax = 20,
                                probRegion=0.95, legend=TRUE)



plot_values  <- matrix(0, ncol=2, nrow=length(obs_to_plot))
count  <- 1
for(k in 1:length(obs_to_plot)) {
  cov_true  <- V  %*% (cov_list[[obs_to_plot[[k]]]]) %*% t(V) +
    error_sd^2  * (Vnull %*% t(Vnull))
  pmPsi  <- t(Vstar)  %*% cov_true  %*% Vstar
  eigenPsi  <- eigen(pmPsi)
  plot_values[count, 1]  <- eigenPsi$values[1]/eigenPsi$values[2]
  plot_values[count, 2]  <- atan(eigenPsi$vectors[2, 1]/eigenPsi$vectors[1, 1])
  count <- count + 1
}

colnames(plot_values)  <- c("y", "x")


## cols <- rev(colorspace::sequential_hcl(4, "Viridis"))

posterior_plot + geom_point(aes(x=x, y=log2(y), col=as.factor(1:3)), shape=18, data=as_tibble(plot_values), size=4)


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


Sig_x_lst  <- lapply(1:n, function(j) {
    j_lst  <- lapply(1:1000, function(i) {
        gx  <- estep$samples$gamma[, , i] * X[j, ]
        gx %*% t(gx)
    })
    Sig_x <- Reduce('+', j_lst)/1000
    Sig_x
})
plot(apply(abs(estep$samples$gamma), c(2, 3), mean), abs(gamma))

