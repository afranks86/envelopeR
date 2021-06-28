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

targeted  <- read_csv("data/untargeted.csv")
subject_info  <- read_csv("data/subject_info.csv")

Y  <- targeted %>% as.matrix

## Remove extreme outliers (more than 4 median asbolute deviations from 0)
Y  <- apply(Y, 2, function(x) {
  mad_val <- mad(x, na.rm=TRUE)
  outliers  <- which(abs(x) > 4*mad_val)
  x[outliers]  <- NA
  x
})
Y[Y==0] <- NA

## Remove metaboltes with a lot of missingess
indices  <- which(colMeans(is.na(Y)) < 0.05)
Y  <- Y[, indices]

## Impute missing values on remaining
Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
Y <- t(Yt) %>% as.matrix

subject_info  <- subject_info %>% mutate(Type2 = ifelse(Type != "AD", "C", "AD"))
X <- subject_info %>% mutate(Type2 = ifelse(Type != "AD", "C", "AD")) %>%
  model_matrix(~ Age + Type2 + Sex) %>%
  as.matrix

## Focus only on aging among controls
control_indices  <- which(subject_info$Type %in% c("CY", "CM", "CO"))

Xfit  <- X[control_indices, c("Age", "SexM"), drop=FALSE]
Yfit  <- Y[control_indices, ]

Xfit[order(Xfit[, 1]), ]
xord  <- order(Xfit[, 1])

## Get the rank
s  <- getRank(Yfit)
q <- ncol(Xfit)

#############################
#### Envelope Fit
############################

res <- fit_envelope(Y=Yfit, X=Xfit, s=s, distn="covreg", maxIters=1000, Vinit="OLS")

YVfit  <- Yfit %*% res$V

covreg_fit  <- covreg::covreg.mcmc(YVfit ~ Xfit - 1,
                                   YVfit ~ Xfit,
                                   R=5, niter=10000,
                                   nthin=10)

cov_psamp  <- covreg::cov.psamp(covreg_fit)
a_psamp  <- covreg_fit$A.psamp

## Pathway analysis for age and sex
for(k in 1:2) {
  if(k == 1)
    type  <- "age"
  else
    type  <- "sex"

  mean_coefs  <- res$V %*% covreg_fit$B1.psamp[, k, ]

  met_info  <- subject_data %>% dplyr::select(Metabolite, `m/z`, `Retention time (min)`, Mode) %>%
    distinct
  met_info_pos   <- met_info %>% filter(Mode == "pos")
  met_info_neg   <- met_info %>% filter(Mode == "neg")

  fname  <- sprintf("mean_pathways_neg_%s.tsv", type)
  out_name  <- sprintf("mean_pathways_neg_%s_out", type)

  neg_indices  <- which(colnames(Yfit) %in% met_info$Metabolite[met_info$Mode == "neg"])
  rownames(mean_coefs)  <- colnames(Yfit)
  mean_coefs_neg  <- mean_coefs[neg_indices, ]
  run_mummichog_mean(mean_coefs_neg, met_info_neg, fname, out_name, mode="negative", work_dir = work_dir)

  fname  <- sprintf("mean_pathways_pos_%s.tsv", type)
  out_name  <- sprintf("mean_pathways_pos_%s_out", type)

  pos_indices  <- which(colnames(Yfit) %in% met_info$Metabolite[met_info$Mode == "pos"])
  rownames(mean_coefs)  <- colnames(Yfit)
  mean_coefs_pos  <- mean_coefs[pos_indices, ]
  run_mummichog_mean(mean_coefs_pos, met_info_pos, fname, out_name, mode="positive", work_dir = work_dir)


  wd_files  <- dir(work_dir)
  sub_dir_pos  <- wd_files[grep(sprintf("mean_pathways_pos_%s_out", type), wd_files)][1]

  pathways_pos  <- read_tsv(paste0(work_dir, sub_dir_pos,
                                   sprintf("/tables/mcg_pathwayanalysis_mean_pathways_pos_%s_out.tsv", type))) %>%
    filter(`p-value` < 0.05, overlap_size > 0) %>% slice(1:5)

  sub_dir_neg  <- wd_files[grep(sprintf("mean_pathways_neg_%s_out", type), wd_files)][1]
  pathways_neg  <- read_tsv(paste0(work_dir, sub_dir_neg,
                                   sprintf("/tables/mcg_pathwayanalysis_mean_pathways_neg_%s_out.tsv", type))) %>%
    filter(`p-value` < 0.05, overlap_size > 0) %>% slice(1:5)

  pathways  <- bind_rows(pathways_pos, pathways_neg)

  tibble(pathway=pathways$pathway, pval = pathways$`p-value`) %>%
    arrange(pval) %>% kable(format="latex") %>% cat(., file=sprintf("untargeted_%s_mean.tex", type))
}

## aging_to_plot
plot_type  <- "female"
X_unique <- unique(Xfit[, c("Age", "GenderM")])

if(plot_type == "female")
  to_plot  <- c(42, 46, 5, 59, 54, 27)
else
  to_plot  <- c(13, 21, 63, 28, 1, 52)
X_unique[to_plot, ]

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")
nms <- X_unique[, 1]


rotV <- rotate_basis(res$V, cov_psamp, to_plot[1], to_plot[length(to_plot)])$rotV
rownames(rotV)  <- colnames(Y)

save(Yfit, Xfit, covreg_fit, cov_psamp, res, file = paste0("untargeted_covreg-", today(), ".Rdata"))

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))

run_mummichog_script  <- TRUE
for(i in seq(1, s, 2)) {

  work_dir  <- "../results/untargeted_pathways/"

  labels  <- pos  <- xy  <- c()
  for(j in i:(i+1)) {
    if(run_mummichog_script) {
      fname  <- sprintf("cov_pathways%i_%s_%s.tsv", j, "neg", plot_type)
      out_name  <- sprintf("cov_pathways_out%i_%s_%s", j, "neg", plot_type)
      rotVneg  <- rotV[intersect(rownames(rotV),
                                    met_info$Metabolite[met_info$Mode == "neg"]), ]
      run_mummichog(rotVneg[, j], met_info_neg, fname, out_name, mode = "negative", work_dir=work_dir, cutoff = 0.05)

      fname  <- sprintf("cov_pathways%i_%s_%s.tsv", j, "pos", plot_type)
      out_name  <- sprintf("cov_pathways_out%i_%s_%s", j, "pos", plot_type)
      rotVpos  <- rotV[intersect(rownames(rotV),
                                 met_info$Metabolite[met_info$Mode == "pos"]), ]
      run_mummichog(rotVpos[, j], met_info_pos, fname, out_name, mode="positive", work_dir=work_dir, cutoff = 0.05)


    }
    wd_files  <- dir(work_dir)
    sub_dir_neg  <- wd_files[grep(sprintf("cov_pathways_out%i_neg_%s", j, plot_type), wd_files)][1]
    pathways_neg  <- read_tsv(paste0(work_dir, sub_dir_neg,
                                     sprintf("/tables/mcg_pathwayanalysis_cov_pathways_out%i_%s_%s.tsv", j, "neg", plot_type))) %>%
    filter(`p-value` < 0.05, overlap_size > 0) %>% slice(1:5)
    sub_dir_pos  <- wd_files[grep(sprintf("cov_pathways_out%i_pos_%s", j, plot_type), wd_files)][1]
    pathways_pos  <- read_tsv(paste0(work_dir, sub_dir_pos,
                                     sprintf("/tables/mcg_pathwayanalysis_cov_pathways_out%i_%s_%s.tsv", j, "pos", plot_type))) %>%
      filter(`p-value` < 0.05, overlap_size > 0) %>% slice(1:5)


    pathways  <- bind_rows(pathways_neg, pathways_pos)
    labs  <- pathways %>% .$pathway %>% unique
    labels  <- c(labels, labs)
    pos <- c(pos, rep(max(rotV[, j]), length(labs)))
    xy  <- c(xy, ifelse(j==rep(i, length(labs)), "x", "y"))
  }


  custom_label_data <- tibble(x=c(pos[xy=="x"], rep(0, sum(xy=="y"))),
                           y=c(rep(0, sum(xy=="x")), pos[xy=="y"]),
                           labels = labels)


  post  <- create_plots(res$V, cov_psamp,
                        n1=to_plot[1], n2=to_plot[length(to_plot)],
                        to_plot = to_plot, obs_names=Xfit[to_plot, "Age"], col_values=cols,
                        labels=colnames(Y), plot_type="posterior", alpha=0.5)

  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols,
                         obs_names=Xfit[to_plot, "Age"], view=c(i, i+1), labels=1:nrow(res$V),
                         legend.pos="right", main="Untargeted",
                         custom_label_data = custom_label_data, label_size=3)
  combo
                                        #colnames(Y))
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_biplot-%i%i.pdf", i, i+1), bp)
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_posterior-%i%i.pdf", i, i+1), post)
  ggsave(sprintf("../figs/untargeted/aging_envelopeR/aging_reg_combo-%i%i_%s.pdf", i, i+1, plot_type), combo, width=14)
}


## sex_to_plot

## Difference in covariance for M/F for 57 year old
plot_type  <- "sex_57"
X_unique <- unique(Xfit[, c("Age", "GenderM")])

to_plot  <- c(43, 32)
X_unique[to_plot, ]

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")
nms <- X_unique[, 2]

rotV <- rotate_basis(res$V, cov_psamp, to_plot[1], to_plot[length(to_plot)])$rotV

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))

run_mummichog_script  <- TRUE
for(i in seq(1, s, 2)) {

  work_dir  <- "../results/untargeted_pathways/"

  labels  <- pos  <- xy  <- c()
  for(j in i:(i+1)) {
    if(run_mummichog_script) {
      fname  <- sprintf("cov_sex_pathways%i_%s_%s.tsv", j, "neg", plot_type)
      out_name  <- sprintf("cov_sex_pathways_out%i_%s_%s", j, "neg", plot_type)
      rotVneg  <- rotV[intersect(rownames(rotV),
                                    met_info$Metabolite[met_info$Mode == "neg"]), ]
      run_mummichog(rotVneg[, j], met_info_neg, fname, out_name, mode = "negative", work_dir=work_dir, cutoff = 0.05)

      fname  <- sprintf("cov_sex_pathways%i_%s_%s.tsv", j, "pos", plot_type)
      out_name  <- sprintf("cov_sex_pathways_out%i_%s_%s", j, "pos", plot_type)
      rotVpos  <- rotV[intersect(rownames(rotV),
                                 met_info$Metabolite[met_info$Mode == "pos"]), ]
      run_mummichog(rotVpos[, j], met_info_pos, fname, out_name, mode="positive", work_dir=work_dir, cutoff = 0.05)


    }
    wd_files  <- dir(work_dir)
    sub_dir_neg  <- wd_files[grep(sprintf("cov_sex_pathways_out%i_neg_%s", j, plot_type), wd_files)][1]
    pathways_neg  <- read_tsv(paste0(work_dir, sub_dir_neg,
                                     sprintf("/tables/mcg_pathwayanalysis_cov_sex_pathways_out%i_%s_%s.tsv", j, "neg", plot_type))) %>%
    filter(`p-value` < 0.05, overlap_size > 0) %>% slice(1:5)
    sub_dir_pos  <- wd_files[grep(sprintf("cov_sex_pathways_out%i_pos_%s", j, plot_type), wd_files)][1]
    pathways_pos  <- read_tsv(paste0(work_dir, sub_dir_pos,
                                     sprintf("/tables/mcg_pathwayanalysis_cov_sex_pathways_out%i_%s_%s.tsv", j, "pos", plot_type))) %>%
      filter(`p-value` < 0.05, overlap_size > 0) %>% slice(1:5)

    if(nrow(pathways_neg == 0) & nrow(pathways_pos == 0))
      next
    else if(nrow(pathways_neg == 0))
      pathways  <- pathways_pos
    else if(nrow(pathways_pos == 0))
      pathways  <- pathways_neg
    else
      pathways  <- bind_rows(pathways_neg, pathways_pos)
    labs  <- pathways %>% .$pathway %>% unique
    labels  <- c(labels, labs)
    pos <- c(pos, rep(max(rotV[, j]), length(labs)))
    xy  <- c(xy, ifelse(j==rep(i, length(labs)), "x", "y"))
  }


  custom_label_data=tibble(x=c(pos[xy=="x"], rep(0, sum(xy=="y"))),
                           y=c(rep(0, sum(xy=="x")), pos[xy=="y"]),
                           labels = labels)

  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols,
                         obs_names=nms[to_plot], view=c(i, i+1), labels=1:nrow(res$V),
                         legend.pos="right", main="Untargeted",
                         custom_label_data = custom_label_data, label_size=3)
                                        #colnames(Y))
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_biplot-%i%i.pdf", i, i+1), bp)
  ## ggsave(sprintf("../figs/gotms/aging_mgCov/aging_reg_posterior-%i%i.pdf", i, i+1), post)
  ggsave(sprintf("../figs/untargeted/aging_envelopeR/sex_reg_combo-%i%i_%s.pdf", i, i+1, plot_type), combo, width=14)
}
