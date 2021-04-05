## envelope model for CSF data
##
library(tidyverse)
library(ggExtra)
library(mvtnorm)
library(Amelia)
library(ggridges)
library(envelopeR)
library(modelr)
library(Matrix)
library(MetaboAnalystR)
library(kableExtra)

source("utility.R")

## load("preprocessed_gotms_data.RData")
##load("preprocessed_gotms_data-2018-03-21.RData")

load("../data/preprocessed_csf_data.RData")

subject_data$Metabolite %>% unique


## X <- matrix(rnorm(N, sd = 5), ncol=1)
wide_data <- subject_data %>%
  filter(!(Type %in% c("Other"))) %>%
  mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Code", "Mode", "RunIndex",
                          "Raw", "Trend", "RawScaled")) %>%
  spread(key = Metabolite, value = Abundance)

Y <- wide_data %>%
    dplyr::select(-one_of("Id", "Type", "Type2", "Gender",
                          "Age", "APOE", "Batch", "Index",
                          "GBAStatus", "GBA_T369M", "cognitive_status")) %>%
    as.matrix()

## Remove extreme outliers (more than 4 median asbolute deviations from 0)
Y  <- apply(Y, 2, function(x) {
  mad_val <- mad(x, na.rm=TRUE)
  outliers  <- which(abs(x) > 4*mad_val)
  x[outliers]  <- NA
  x
})


## Impute missing values using
Y <- amelia(Y, m = 1, empri = 100)$imputations$imp1


X <- wide_data %>%
  dplyr::select(one_of("Gender", "Age", "Type2", "APOE")) %>%
  mutate(APOE = factor(APOE, c("22", "23", "24", "33", "34", "44"))) %>%
  model_matrix(~ Age + Type2 + Gender + APOE) %>% as.matrix

## Focus only on aging among controls
control_indices  <- which(wide_data$Type2=="C")

Xfit  <- X[control_indices, c("Age", "GenderM"), drop=FALSE]
Yfit  <- Y[control_indices, ]

Xfit[order(Xfit[, 1]), ]
xord  <- order(Xfit[, 1])

glm_res  <- glmnet(Xfit, Yfit, family="mgaussian")

D <- diag(1/sqrt(apply(Y, 2, function(x) mad(x, na.rm=TRUE))))
# D <- diag(ncol(Yfit))

s <- getRank(Yfit)
r <- 0
q <- ncol(X)

fit_obj  <- lm(Yfit ~ Xfit)
getRank(fit_obj$residuals)

indices  <- 1:nrow(Xfit)
prior_counts  <- 1

res <- fit_envelope(Y=Yfit, X=scale(Xfit[indices, ]), D=D, s=s, r=0, distn="covreg",
                               prior_counts=prior_counts,
                               maxIters=1000, U1=0*diag(s), U0=0*diag(r),
                               Vinit="OLS", L=0)

covreg_fit  <- covreg::covreg.mcmc(YVfit ~ scale(Xfit[indices, ]) - 1,
                                   YVfit ~ scale(Xfit[indices, ]),
                                   R=10, niter=10000, nthin=10)

## mean fit
mean_coefs  <- covreg_fit$B1.psamp

## Covariance Fit
cov_psamp  <- covreg::cov.psamp(covreg_fit)
mean_coefs_age  <- res$V %*% mean_coefs[, 1, ]
mean_coefs_sex  <- res$V %*% mean_coefs[, 2, ]
rownames(mean_coefs_age)  <- rownames(mean_coefs_sex)  <-  colnames(Y)

map_dfr(list(Age=mean_coefs_age, Sex=mean_coefs_sex), function(mat) apply(mat, 1,
                     function(x) {
                       frac_neg  <- mean(x < 0)
                       pval <- 2*min(frac_neg, 1-frac_neg)
                       tstat  <- mean(x) / sd(x)
                       c("P-value"=pval, "T-statistic"=tstat)
                     }) %>% t %>% as_tibble(rownames="Metabolite"), .id="Type") %>%
  group_by(Type) %>% 
  arrange(`P-value`, desc(abs(`T-statistic`))) %>%
  mutate(`Q-value` = `P-value`*n()/row_number()) %>%
  ungroup() ->
  regression_stats


regression_stats %>% filter(`Q-value` < 0.05, Type == "Sex") %>%
  kable(format="latex") %>%
  cat(., file = "targeted_sex.tex")

regression_stats %>% filter(`Q-value` < 0.05, Type == "Age") %>%
  kable(format="latex") %>%
  cat(., file = "targeted_age.tex")

index_map  <- match(apply(Xfit, 1, function(x) paste(x, collapse="_")), apply(X_unique, 1, function(x) paste(x, collapse="_")))

compute_log_likelihood(YVfit, index_map,
                       apply(covreg::m.psamp(covreg_fit), 1:2, mean),
                       apply(cov_psamp, 1:3, mean))


X_unique <- unique(Xfit[indices, c("Age", "GenderM")])

plot(X_unique[, "Age"], rowMeans(apply(cov_psamp, c(1, 4), function(x) norm(cov2cor(x), type="F"))), col=(X_unique[, "GenderM"]+1))

rowMeans(apply(cov_psamp, c(1, 4), function(x) norm(cov2cor(x), type="F")))
apply(apply(cov_psamp, c(1, 4), function(x) norm(cov2cor(x), type="F")), 1, sd)


as_tibble(apply(cov_psamp, c(1, 4), function(x) norm(cov2cor(x), type="F"))) %>%
  mutate(age = X_unique[, "Age"], sex=X_unique[, 2])  %>%
  pivot_longer(cols = -c("age", "sex") , names_to="sample", values_to="norm") %>%
  filter(sex==1) %>%
  ggplot() + geom_line(aes(x=age, y=norm, group=sample), alpha=0.2) + theme_bw()


X_unique <- unique(Xfit[, c("Age", "GenderM")])

cbind(X_unique, 1:nrow(X_unique))[order(X_unique[, 1]), ]

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")
nms <- X_unique[, 1]



## aging_to_plot
to_plot  <- c(13, 21, 63, 28, 1, 52) ## MALES
to_plot  <- c(44, 26, 43, 2, 27) ## FEMALES
X_unique[to_plot, ]

save(Yfit, Xfit, covreg_fit, cov_psamp, Vfit, res, file = paste0("targeted_covreg-", today(), ".Rdata"))


cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))
for(i in seq(1, 29, 2)) {

  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols, nlabeled=20,
                         obs_names=nms, view=c(i, i+1), labels=colnames(Y), main="Targeted", legend.title="Age")
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_biplot-%i%i.pdf", i, i+1), bp)
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_posterior-%i%i.pdf", i, i+1), post)
  ggsave(sprintf("../figs/targeted/aging_envelopeR/aging_reg_combo-%i%i.pdf", i, i+1), combo, width=14)
}

## gender_to_plot
to_plot  <- c(43, 32)
X_unique[to_plot, ]
nms <- X_unique[, 2]

cbind(X_unique, 1:nrow(X_unique))[order(X_unique[, 1]), ]
to_plot  <- c(49, 14, 15, 58)
to_plot  <- c(14, 49, 58, 15)
X_unique[to_plot, ]

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))
for(i in seq(1, s, 2)) {

  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols,
                         obs_names=nms, view=c(i, i+1), labels=colnames(Y))
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_biplot-%i%i.pdf", i, i+1), bp)
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_posterior-%i%i.pdf", i, i+1), post)
  ggsave(sprintf("../figs/targeted/aging_envelopeR/sex_reg_combo-%i%i.pdf", i, i+1), combo, width=14)
}

## gender_to_plot
cbind(X_unique, 1:nrow(X_unique))[order(X_unique[, 1]), ]
to_plot  <- c(49, 14, 15, 58)
to_plot  <- c(14, 49, 58, 15)
X_unique[to_plot, ]

for(i in seq(1, 11, 2)) {
  ## bp  <- create_plots(res$V, cov_psamp, n1=to_plot[1], n2=to_plot[length(to_plot)], to_plot = to_plot,
  ##                     obs_names=nms, view=c(i, i+1), plot_type="biplot")
  ## post  <- create_plots(Vaging, cov_psamp, n1=to_plot[1], n2=to_plot[length(to_plot)], to_plot = to_plot,
  ##                       obs_names=nms, view=c(i, i+1), plot_type="posterior")

  combo  <- create_plots(res$V, cov_psamp, n1=to_plot[1], n2=to_plot[length(to_plot)], to_plot = to_plot,
                         obs_names=nms, view=c(i, i+1), labels=colnames(Y))
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_biplot-%i%i.pdf", i, i+1), bp)
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_posterior-%i%i.pdf", i, i+1), post)
  ggsave(sprintf("../figs/targeted/aging_envelopeR/gender_reg_combo-%i%i.pdf", i, i+1), combo, width=14)
}



## Save all metabolites
all_mets <- subject_data$Metabolite %>% unique %>% as.data.frame
write_csv(all_mets, "all_mets.csv")

i  <- 7
rotation_result <- rotate_basis(res$V, cov_psamp, n1=to_plot[1], n2=to_plot[length(to_plot)])
V_rotated  <- rotation_result$rotV[, i:(i+1)]
rownames(V_rotated) <- colnames(Y)
sig_mets1 <- names(sort(abs(V_rotated[, 1]), decreasing=TRUE))[1:20]
sig_mets2 <- names(sort(abs(V_rotated[, 2]), decreasing=TRUE))[1:20]
sig_mets_both <- names(sort(rowSums(V_rotated[, 1:2]^2), decreasing=TRUE))[1:20]

Vall <- rowSums(rotation_result$rotV^2)
names(Vall) <- colnames(Y)
sig_mets_both <- names(sort(Vall, decreasing=TRUE)[1:10])

mSet <- InitDataObjects("conc", "msetora", FALSE)
mSet <- Setup.HMDBReferenceMetabolome(mSet, "all_mets.csv");
cmpd.vec <- sig_mets_both
mSet <- Setup.MapData(mSet, cmpd.vec);
mSet <- CrossReferencing(mSet, "name");
mSet <- CreateMappingResultTable(mSet)
mSet <- SetMetabolomeFilter(mSet, T);
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet <- CalculateHyperScore(mSet)
mSet <- CalculateHyperScore(mSet)

## Scores
mSet$analSet$ora.mat
mSet$name.map$map.table
mSet$analSet$ora.hits

as_tibble(mSet$dataSet$map.table)$KEGG

