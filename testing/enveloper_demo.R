library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(mvtnorm)
library(covreg)
library(envelopeR)


test_data <- generate_test_data(q=2, cov_rank=2, intercept=TRUE, gamma_sd=1, error_sd=1, seed=4445)
Y <- test_data$Y
X  <- test_data$X
s  <- test_data$s
p  <- test_data$p
V  <- test_data$V

envfit <- fit_envelope(Y, X, distn="covreg", s=s,
                       Vinit="OLS",
                       verbose_covreg=FALSE)

tr(envfit$V  %*% t(envfit$V)  %*% V  %*% t(V))/s

YV  <- Y %*% envfit$V

#cov_psamp  <- covreg::cov.psamp(envfit$covariance_list$covreg_res)

res  <-  covreg::covreg.mcmc(YV ~ X - 1, YV ~ X, niter=10000, nthin=10)
cov_psamp  <- covreg::cov.psamp(res)
cov_psamp  <- cov_psamp[, , , 1:1000]


type  <- "mag"
## obs_to_plot <- sample(100, 3)
ix <- sort(X[, 1], index.return=TRUE)$ix
obs_to_plot <- ix[c(1, 25, 50, 75, 100)]
names(obs_to_plot)  <- c(0, 0.25, 0.5, 0.75, 1)

cols <- colorspace::sequential_hcl(5, "viridis")


post  <- create_plots(envfit$V, cov_psamp,
                      n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                      to_plot = obs_to_plot, col_values=cols,
                      labels=colnames(Y), plot_type="posterior")
post
