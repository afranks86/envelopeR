This folder includes scripts to generate all of the simulation plots included in Biometrics paper.

All figures can be generated using the `plot_simulation_results.R` script. All of the results are generated using the following scripts:

- Results for Figure 1a can be generated from `subspace_misspecification_1a.R`.  The generated plot shows the percentage change in risk for envelopeR estimates using the wrong subspace dimension.
- Results for Figure 1b can be generated from `covariates_1b.R`.  The generated plot illustrates the percent increase in risk of the naive estimator over the envelopeR estimator as a function of the number of predictors.  
- Results for Figure 2 can be generated from `goodness_of_fit_2.R`.  This shows the goodness of fit of the models with subspace dimension which is smaller than the true subspace dimension.  
- Results for Online supplimentary Figure 1 can be generated from `initialization_S1.R`.  Demonstrates the importance of carefully chosen initial subspaces (random initializations find suboptimal solutions).
- Results Online supplementary Figure 2 can be generated from `scalability_S2.R` .  We recommend using the cloud, as the number of replicates, features and sample sizes fit in this simulation mean that the script takes several days to complete.
