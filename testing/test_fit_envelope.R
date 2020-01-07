

library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
source("utility_functions.R")
source("envelope_functions.R")

n <- 200
p <- 800
s <- 40
r <- 0
q <- 5

## n <- 1000
## p <- 500
## s <- 50
## r <- 10
## q <- 10

X <- matrix(rnorm(n*q), nrow=n, ncol=q)

V <- rustiefel(p, s)
V <- rbind(diag(s), matrix(0, nrow=p-s, ncol=s))

eta <- matrix(rnorm(q*s, sd=0.1), nrow=q)
beta <- eta %*% t(V)

YV <- X %*% eta + rmvnorm(n, mean= rep(0, s),
                          sigma = diag(sort(rexp(s, 1/8), decreasing=TRUE)))

if(r > 0) {
  U <- NullC(V) %*% rustiefel(p-s, r)
  YU <- rmvnorm(n, mean= rep(0, r),
                sigma = diag(sort(rexp(r, 1/8), decreasing=TRUE)))
} else {
  U <- matrix(nrow=p, ncol=0)
  YU <- matrix(0, nrow=n, ncol=0)
}

Y <- YV %*% t(V) + YU %*% t(U) +
  matrix(rnorm(n * p, sd=3), nrow=n, ncol=p)


res <- fit_envelope(Y, X, s=s, r=0,
                    distn="normal",
                    nchunks=1,
                    Vinit="COV", use_py=FALSE,
                    maxIters=10000,
                    searchParams=list(rho=1e-5, eta=0.9),
                    L = 50,
                    prior_counts = 50)
res$F(cbind(V, U))
res$F(rustiefel(p, s+r))
res$F(res$V)
tr(res$V %*% t(res$V) %*% cbind(V, U) %*% t(cbind(V, U)))  / (s+r)

sum((res$beta_ols - beta)^2)
sum((res$beta_env - beta)^2)
mean(rowSums(res$V^2) < 1e-10)
mean(rowSums(cbind(V, U)) < 1e-10)

#######################################################################

## Covariance Regression model

## Generate some fake data
library(mvtnorm)
P <- 50
N <- 50
S <- 20
R <- 20
Q <- 2

## Generate envelope
V <- rustiefel(P, S)
Vperp <- NullC(V)

## regression parameters

## eta small relative to evals gives based improvements in MSE

## eta <- matrix(c(-1, 1), ncol=1)

eta <- matrix(c(-10, 10, 10, 2), nrow=2, ncol=2)
eta <- matrix(rnorm(S^2, 0, 1), nrow=S)
Beta <- t(eta) %*% t(V)




X <- matrix(rnorm(N*Q), nrow=N, ncol=Q)

eta <- matrix(rnorm(Q*S, sd=2), nrow=Q)

beta <- matrix(rnorm(Q*S, sd=2), nrow=Q)


get_cov <- function(X, beta) {
  beta_X <- X %*% beta
  sig_x <- t(beta_X) %*%  beta_X
  sig_x
}


SigInvXList <- lapply(1:N, function(i) solve(get_cov(X[i, ], beta) + diag(S)))
etaSigXList <- lapply(1:N, function(i) eta %*% SigInvXList[[i]])

Y <- sapply(1:N, function(i)
  rmvnorm(1,
          mean = X[i, ] %*% eta %*% t(V),
          sigma = V %*% solve(SigInvXList[[i]]) %*% t(V) + diag(P))) %>% t
svd(Y)$d


## covariates

X <- matrix(rmvnorm(N, rep(0, S), sigma=diag(rep(5^2, S))), ncol=S)
                                        #X <- rmvnorm(N, mean=c(0,0), sigma=diag(c(5,5)))

## covariance matrices
psi <- diag(sort(rexp(S, 1 / 100), decreasing = TRUE))
## psi <- diag(c(10, 10))
## psi.0 <- diag(rep(100, P-S))
psi.0 <- diag(sort(rexp(P - S, 1 / 100), decreasing = TRUE))

ecdf(diag(psi.0))(diag(psi))

Sigma <- V %*% psi %*% t(V) + Vperp %*% psi.0 %*% t(Vperp)

Y <- sapply(1:N, function(i) rmvn(1, X[i, ] %*% Beta, Sigma)) %>% t

betaHat <- solve(t(X) %*% X) %*% t(X) %*% Y

R <- 2
betav <- svd(betaHat)$v
Vinit <- cbind(betav, NullC(betav)[, 1:R])
prior_counts <- 10




X <- matrix(rnorm(n*q), nrow=n, ncol=q)

V <- rustiefel(p, s)
V <- rbind(diag(s), matrix(0, nrow=p-s, ncol=s))

eta <- matrix(rnorm(q*s, sd=0.1), nrow=q)
beta <- eta %*% t(V)

YV <- X %*% eta + rmvnorm(n, mean= rep(0, s),
                          sigma = diag(sort(rexp(s, 1/8), decreasing=TRUE)))

if(r > 0) {
  U <- NullC(V) %*% rustiefel(p-s, r)
  YU <- rmvnorm(n, mean= rep(0, r),
                sigma = diag(sort(rexp(r, 1/8), decreasing=TRUE)))
} else {
  U <- matrix(nrow=p, ncol=0)
  YU <- matrix(0, nrow=n, ncol=0)
}

Y <- YV %*% t(V) + YU %*% t(U) +
  matrix(rnorm(n * p, sd=3), nrow=n, ncol=p)


res <- fit_envelope(Y, X, s=s, r=0,
                    distn="normal",
                    nchunks=1,
                    Vinit="COV", use_py=FALSE,
                    maxIters=10000,
                    searchParams=list(rho=1e-5, eta=0.9),
                    L = 50,
                    prior_counts = 50)
res$F(cbind(V, U))
res$F(rustiefel(p, s+r))
res$F(res$V)
tr(res$V %*% t(res$V) %*% cbind(V, U) %*% t(cbind(V, U)))  / (s+r)

sum((res$beta_ols - beta)^2)
sum((res$beta_env - beta)^2)
mean(rowSums(res$V^2) < 1e-10)
mean(rowSums(cbind(V, U)) < 1e-10)


#################################################################


res_kd <- fit_envelope(Y, X, s=s, r=r,
                       distn="normal",
                       nchunks=3,
                       prior_counts = 0,
                       Vinit="OLS", use_py=FALSE,
                       maxIters=1000,
                       searchParams=list(rho=1e-5, eta=0.9))

sum((res_kd$beta_ols - beta)^2)
sum((res_kd$beta_env - beta)^2)


res_cv <- cv.glmnet(X, Y, family="mgaussian", alpha=0, intercept=FALSE)
res <- glmnet(X, Y, family="mgaussian", lambda=res_cv$lambda.min, alpha=0, intercept=FALSE)
sum((do.call(cbind, coef(res))[-1, ] - beta)^2)


res_kd$F(cbind(V, U))
res_kd$F(rustiefel(p, s+r))
res_kd$F(res_kd$V)
res_kd$F(res$V)
res$F(cbind(V, U))


vbeta <- svd(beta)$v

tr(res_kd$V[, 1:s] %*% t(res_kd$V[, 1:s]) %*% V %*% t(V)) / s
tr(res_kd$V[, 1:s] %*% t(res_kd$V[, 1:s]) %*% vbeta %*% t(vbeta)) / ncol(vbeta)
tr(res$V[, 1:s] %*% t(res$V[, 1:s]) %*% V %*% t(V)) / s
tr(res$V[, 1:s] %*% t(res$V[, 1:s]) %*% vbeta %*% t(vbeta)) / ncol(vbeta)

tr(res_kd$V[, (s+1):(s+r)] %*% t(res_kd$V[, (s+1):(s+r)]) %*%
   U %*% t(U)) / r
tr(res_kd$V %*% t(res_kd$V) %*% cbind(V, U) %*% t(cbind(V, U))) / (s+r)

res_kd$F(res_kd$V)
res_kd$F(res$V)



res <- fit_envelope(Y, X, s=s, r=r,
                    distn="normal",
                    Vinit="OLS", use_py=FALSE,
                    maxIters=10000,
                    searchParams=list(rho=1e-4, eta=0.9))



Yproj <- Y %*% res$V
eta_hat_env <- solve(t(X) %*% X) %*% t(X) %*% Yproj
beta_env <- eta_hat_env %*% t(res_kd$V)
sum((beta - beta_env)^2)



res_kd$V

sum((res_kd$beta_ols - beta)^2)
sum((res_kd$beta_env - beta)^2)

res_kd$F(res_kd$V)
res_kd$F(res$V)

res$F(res_kd$V)

res <- fit_envelope(Y, X, s=s, r=r,
                    distn="cook",
                    Vinit="COV", use_py=FALSE,
                    maxIters=10000,
                    searchParams=list(rho=1e-4, eta=0.9))

sum((res$beta_ols - beta)^2)
sum((res$beta_env - beta)^2)

res_kd <- fit_envelope(Y, X, s=s, r=r,
                       distn="cook_kd",
                       Vinit="OLS", use_py=FALSE,
                       maxIters=1000,
                       searchParams=list(rho=1e-4, eta=0.9))
Vhat1 <- res_kd$V[, 1:s, drop=FALSE]
res_kd$F(res_kd$V)
res_kd$F(V)

Yproj <- Y %*% res$V
eta_hat_env <- solve(t(X) %*% X) %*% t(X) %*% Yproj
beta_env <- eta_hat_env %*% t(res$V)
sum((beta - beta_env)^2)


tr(t(V) %*% Vhat1 %*% (t(Vhat1) %*% V)) / s
tr(t(U) %*% Vhat1 %*% (t(Vhat1) %*% U))
tr(t(V) %*% Vhat2 %*% (t(Vhat2) %*% V)) / 2
tr(t(V) %*% Vhat %*% (t(Vhat) %*% V)) / 4

Vinit <- cbind(res_kd$V, NullC(res_kd$V)[, 1:r])
res_kd <- fit_envelope(Y, X, s=50, r=10,
                       distn="cook_kd",
                       Vinit=Vinit, use_py=FALSE,
                       maxIters=1000,
                       searchParams=list(rho=1e-4, eta=0.9))

res_kd$F(res_kd$V)
res_kd$F(V)





res <- fit_envelope(Y, X, s=1, r=1,
                    Vinit="OLS", use_py=FALSE,
                    maxIters=1000,
                    searchParams=list(rho=1e-4, eta=0.9))



Vhat1 <- res$V[, 1:s, drop=FALSE]
NV <- NullC(Vhat1)

res2 <- fit_envelope(Y %*% (diag(p) - Vhat1 %*% t(Vhat1)), X, s=1, r=19,
                     distn="cook",
                     Vinit="OLS", use_py=FALSE,
                     maxIters=1000,
                     searchParams=list(rho=1e-4, eta=0.9))
Vhat2 <- res2$V[, 1:2]

t(Vhat1) %*% Vhat2
Vhat <- cbind(Vhat1, Vhat2)
tr(t(V) %*% Vhat1 %*% (t(Vhat1) %*% V)) / s
tr(t(U) %*% Vhat1 %*% (t(Vhat1) %*% U))
tr(t(V) %*% Vhat2 %*% (t(Vhat2) %*% V)) / 2
tr(t(V) %*% Vhat %*% (t(Vhat) %*% V)) / 4


beta_hat <- solve(t(X) %*% X) %*% (t(X) %*% Y)
resid <- Y - X %*% beta_hat
res_proj <- (diag(p) - Vhat %*% t(Vhat)) %*% t(resid)
Uhat <- svd(res_proj)$u[, 1:r]

res3 <- fit_envelope(Y, X, s=s, r=r,
                     Vinit=cbind(Vhat, Uhat), use_py=FALSE,
                     maxIters=1000,
                     searchParams=list(rho=1e-4, eta=0.9))
res3$F(res3$V)
res3$F(cbind(V, U))




res2 <- fit_envelope(Y, X, s=s, r=r,
                     Vinit=cbind(V, U), use_py=FALSE,
                     maxIters=1000, U1=0.1*diag(s), U0=0.1*diag(r),
                     searchParams=list(rho=0.1, eta=0.9))

names(res)
Vhat <- res$V[, 1:s]
Uhat <- res$V[, (s+1):(s+r)]

res$F(res$V)
res$F(cbind(V, U))
res$F(cbind(V))

tr((t(Vhat) %*% V) %*% (t(V) %*% Vhat)) / s
tr((t(Uhat) %*% U) %*% (t(U) %*% Uhat)) / r

                                        #V contained in res$V but not first s columns thereof

tr((t(V) %*% res$V) %*%  (t(res$V) %*%  V) ) / s
## but not first s columns thereof
tr((t(Vhat) %*% V) %*% (t(V) %*% Vhat)) / s

norm(res$dF(res$V))
norm(res$dF(cbind(V, U)))

resid <- Y - X %*% res$beta_ols

tr((t(U) %*% res$V) %*%  (t(res$V) %*%  U) ) / r
tr((t(cbind(V, U)) %*% res$V) %*%  (t(res$V) %*%  cbind(V, U)) ) / (s+r)


sum((beta - res$beta_ols)^2)
sum((beta - res$beta_env)^2)

Z <- Y %*% res$V
beta_new <- solve(t(X) %*% X) %*% t(X) %*% Z %*% t(res$V)
sum((beta - beta_new)^2)

Z <- Y %*% Vhat
beta_new <- solve(t(X) %*% X) %*% t(X) %*% Z %*% t(Vhat)

sum((beta - beta_new)^2)

Z <- Y %*% V
beta_new <- solve(t(X) %*% X) %*% t(X) %*% Z %*% t(V)
sum((beta - beta_new)^2)



p <- 4000
n <- 20
Sigma = diag(sqrt(p:1))

library(mvtnorm)
Y2 <- rmvnorm(n, sigma = Sigma)

shrinkage_coefs <- lw_shrinkage(Y2)
shrinkage_coefs$b2 / shrinkage_coefs$d2



A1 <- S - Sigma
tr(A1 %*% t(A1))/p

A2 <- Sstar - Sigma
tr(A2 %*% t(A2))/p

n <- 200
p <- 10
s <- 2
r <- 0
q <- 5

## n <- 1000
## p <- 500
## s <- 50
## r <- 10
## q <- 10

X <- matrix(rnorm(n*q), nrow=n, ncol=q)

V <- rustiefel(p, s)
eta <- matrix(rnorm(q*s, sd=2), nrow=q)
beta <- eta %*% t(V)

YV <- X %*% eta + rmvnorm(n, mean= rep(0, s),
                          sigma = diag(sort(rexp(s, 1/4), decreasing=TRUE)))

if(r > 0) {
  U <- NullC(V) %*% rustiefel(p-s, r)
  YU <- rmvnorm(n, mean= rep(0, r),
                sigma = diag(sort(rexp(r, 1/64), decreasing=TRUE)))
} else {
  U <- matrix(nrow=p, ncol=0)
  YU <- matrix(0, nrow=n, ncol=0)
}

Y <- YV %*% t(V) + YU %*% t(U) +
  matrix(rnorm(n * p, sd=2), nrow=n, ncol=p)




envelope_fit <- fit_envelope(Y, X,
                             distn = "normal",
                             s=s, r=0,
                             nchunks=1,
                             Vinit=rustiefel(p, s+r), use_py=FALSE,
                             prior_counts = 0, maxIters=100,
                             search_params=list(rho=1e-5, eta=0.9))









x <- matrix(rnorm(100), ncol=10, nrow=10)
Sig <- t(x) %*% x
g <- rustiefel(10, 2)
g2 <- NullC(g)

a <- det(t(g) %*% Sig %*% g)
b <- det(t(g) %*% solve(Sig) %*% g)
a * det(t(g2) %*% solve(Sig) %*% g2)


log(det(Sig))
log(det(t(g) %*% Sig %*% g)) + log(det(t(g2) %*% Sig %*% g2))
