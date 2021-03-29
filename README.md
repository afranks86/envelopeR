
## Background

`envelopeR` is the R package associated with the paper \`\`Reducing
Subspace Models for Large-Scale Covariance Regression’’.

**Summary:** We develop an envelope model for joint mean and covariance
regression in the large \(p\), small \(n\) setting. In contrast to
existing envelope methods, which improve mean estimates by incorporating
estimates of the covariance structure, we focus on identifying
covariance heterogeneity by incorporating information about mean-level
differences. We use a Monte Carlo EM algorithm to identify a
low-dimensional subspace which explains differences in both means and
covariances as a function of covariates, and then use MCMC to estimate
the posterior uncertainty conditional on the inferred low-dimensional
subspace. We demonstrate the utility of our model on a motivating
application on the metabolomics of aging.

Preprint: <https://arxiv.org/abs/2010.00503>

## Installation

``` r
devtools::install_github("afranks86/envelopeR")

library(tidyverse)
library(covreg)
library(envelopeR)
library(rstiefel)
library(mvtnorm)
```

## A Test Dataset

``` r
test_data <- generate_test_data(q=2, cov_rank=2, intercept=TRUE, gamma_sd=1, error_sd=1, seed=4445)
Y <- test_data$Y
X  <- test_data$X
s  <- test_data$s
p  <- test_data$p
V  <- test_data$V
```

## Fit the Envelope model

``` r
## `get_rank` can be used to infer the appropriate dimension of hte subspace of material variation
s_hat <- getRank(Y)

envfit <- fit_envelope(Y, X, distn="covreg", s=s_hat,
                       Vinit="OLS",
                       verbose_covreg=FALSE)
```

    ## [1] "Starting e-step sampling..."
    ## [1] "Finished e-step sampling..."
    ## [1] "------ F(V) = 26.162593 --------"
    ## [1] "Iteration 1: F = 25.935205, dF = Inf, tau = 0.057668"
    ## [1] "Iteration 2: F = 25.655387, dF = -42.059817, tau = 0.046354"
    ## [1] "Iteration 3: F = 25.510635, dF = -30.753482, tau = 0.048745"
    ## [1] "Iteration 4: F = 25.461795, dF = -25.373999, tau = 0.048504"
    ## [1] "Iteration 5: F = 25.398668, dF = -26.382262, tau = 0.034272"
    ## [1] "Iteration 6: F = 25.380274, dF = -11.474626, tau = 0.024779"
    ## [1] "Iteration 7: F = 25.444786, dF = -3.225030, tau = 0.030877"
    ## [1] "Iteration 8: F = 25.534055, dF = -4.057599, tau = 0.022905"
    ## [1] "Iteration 9: F = 25.553231, dF = -5.833445, tau = 0.003387"
    ## [1] "Iteration 10: F = 25.557681, dF = -6.241172, tau = 0.000741"
    ## [1] "Iteration 11: F = 25.558807, dF = -6.335771, tau = 0.000185"
    ## [1] "Iteration 12: F = 25.559089, dF = -6.359692, tau = 0.000046"
    ## [1] "Iteration 13: F = 25.559094, dF = -6.365691, tau = 0.000001"
    ## [1] "Iteration 14: F = 25.559096, dF = -6.365785, tau = 0.000000"
    ## [1] "Iteration 15: F = 25.559096, dF = -6.365832, tau = 0.000000"
    ## [1] "Iteration 16: F = 25.559096, dF = -6.365837, tau = 0.000000"
    ## [1] "Reached maximum iterations in line search."
    ## [1] "Iteration 17: F = 25.559096, dF = -6.365839, tau = 0.000000"
    ## [1] "F(V) = 25.559096, time = 7.25054812431335"
    ## [1] "Starting e-step sampling..."
    ## [1] "Finished e-step sampling..."
    ## [1] "------ F(V) = 24.626351 --------"
    ## [1] "Iteration 1: F = 24.623553, dF = Inf, tau = 0.064590"
    ## [1] "Iteration 2: F = 24.623589, dF = -0.027768, tau = 0.106081"
    ## [1] "Iteration 3: F = 24.623618, dF = -0.065168, tau = 0.025052"
    ## [1] "Iteration 4: F = 24.623944, dF = -0.015030, tau = 0.027631"
    ## [1] "Iteration 5: F = 24.624090, dF = -0.015763, tau = 0.009292"
    ## [1] "Iteration 6: F = 24.624105, dF = -0.016954, tau = 0.000885"
    ## [1] "Iteration 7: F = 24.624108, dF = -0.017099, tau = 0.000194"
    ## [1] "Iteration 8: F = 24.624109, dF = -0.017132, tau = 0.000096"
    ## [1] "Iteration 9: F = 24.624110, dF = -0.017148, tau = 0.000024"
    ## [1] "Iteration 10: F = 24.624110, dF = -0.017152, tau = 0.000006"
    ## [1] "Iteration 11: F = 24.624110, dF = -0.017153, tau = 0.000000"
    ## [1] "Iteration 12: F = 24.624110, dF = -0.017153, tau = 0.000000"
    ## [1] "Iteration 13: F = 24.624110, dF = -0.017153, tau = 0.000000"
    ## [1] "Iteration 14: F = 24.624110, dF = -0.017153, tau = 0.000000"
    ## [1] "F(V) = 24.624110, time = 4.82364702224731"
    ## [1] "Starting e-step sampling..."
    ## [1] "Finished e-step sampling..."
    ## [1] "------ F(V) = 24.602381 --------"
    ## [1] "Reached maximum iterations in line search."
    ## [1] "Iteration 1: F = 24.602381, dF = Inf, tau = 0.000000"
    ## [1] "F(V) = 24.602381, time = 4.13055491447449"
    ## [1] "Starting e-step sampling..."
    ## [1] "Finished e-step sampling..."

Check inferred subspace is close to the true material subspace (1 means
they are identical subspaces).

``` r
print(sprintf("Subspace similarity: %f", 
              tr(envfit$V  %*% t(envfit$V)  %*% V  %*% t(V))/s))
```

    ## [1] "Subspace similarity: 0.993486"

``` r
## Sphericity test: if small pvalue, reject isotropy of immaterial subspace
Vperp <- NullC(envfit$V) 
print(sprintf("P-value of sphericity test is: %f", 
              test_sphericity(Y %*% Vperp)$p.value))
```

    ## [1] "P-value of sphericity test is: 0.882835"

Re-run full Bayesian inference conditional on the inferred subspace of
material variation. We do this to run the sampler for longer iterations
than the default returned by the EM algorithm. We can also get
`cov_psamp` from the list returned by `fit_envelope` by running
`cov_psamp <- cov.psamp(envfit$covariance_list$covreg_res)` (default 100
iterations).

``` r
YV  <- Y %*% envfit$V
res  <-  covreg::covreg.mcmc(YV ~ X - 1, YV ~ X, niter=10000, nthin=10, verb=FALSE)
cov_psamp  <- covreg::cov.psamp(res)
cov_psamp  <- cov_psamp[, , , 1:1000]
```

## Plotting

### Posterior plot

We plot posterior samples of \(\Psi_x\) for minimum, maximum and
quartiles of the first covariate \(X_1\). Non-overlapping groups is
suggestive of differences in the covariances a posteriori.

``` r
type  <- "mag"

## plot posterior distributions for min, max and quartiles of X
ix <- sort(X[, 1], index.return=TRUE)$ix
obs_to_plot <- ix[c(1, 25, 50, 75, 100)]

## values correspond to the quantiles of x
names(obs_to_plot)  <- c(0, 0.25, 0.5, 0.75, 1)

cols <- colorspace::sequential_hcl(5, "viridis")

post  <- create_plots(envfit$V, cov_psamp,
                      n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                      to_plot = obs_to_plot, col_values=cols,
                      labels=colnames(Y), plot_type="posterior", alpha=0.5)
```

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
    ## Using compatibility `.name_repair`.

``` r
post
```

![](README_files/figure-gfm/posterior_plot-1.png)<!-- -->

The true eigenvalues and angles on this two-dimensional subspaces are:

    ## # A tibble: 5 x 3
    ##   quantile  eval angle
    ##   <chr>    <dbl> <dbl>
    ## 1 0          8.5 1.23 
    ## 2 0.25      10.6 1.00 
    ## 3 0.5       13.3 0.723
    ## 4 0.75      15.3 0.522
    ## 5 1         18.3 0.204

### Biplot

Contours show posterior mean covariances matrices for \(\Psi_x\). Points
represent the largest magnitue loadings on this two-deimsinao subspace.
We choose the subspace to maximize the difference in covariance matrices
between observation \(n_1\) and \(n_2\)

``` r
rownames(envfit$V)  <- 1:p
biplot  <- create_plots(envfit$V, cov_psamp,
                        n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                        to_plot = obs_to_plot, col_values=cols,
                        labels=colnames(Y), plot_type="biplot")
biplot
```

![](README_files/figure-gfm/biplot-1.png)<!-- -->
