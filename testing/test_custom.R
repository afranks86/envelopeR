library(envelopeR)
library(rstan)

sm  <- stan_model("../src/stan_files/mvregression.stan")

sampler <- function(YV, X) {

  stan_fit  <- sampling(sm, data=list(y=YV, x=X, K=ncol(Y), J=ncol(X), N=nrow(YV)))

  SigInvXList  <-  estep$SigInvList
  muSigXList  <-  estep$muSigInvList

}

library(tidyverse)
library(ggridges)
library(rstiefel)
library(modelr)
library(patchwork)
library(envelopeR)
library(cowplot)
library(RColorBrewer)

cattle_data <- read_csv("../testing/cattle.txt")
cattle_data <- cattle_data %>% mutate(Treatment = factor(ifelse(Treatment == "A", "Control", "Treatment")))

cattle_spaghetti <- cattle_data %>% gather(key=Date, value=Weight, -(1:2)) %>% ggplot(aes(x=as.numeric(Date), y=Weight, colour=Treatment, group=id)) + geom_line() + xlab("Days")

cattle_mean <- cattle_data %>% gather(key=Date, value=Weight, -(1:2)) %>% group_by(Treatment, Date) %>% summarise(Mean=mean(Weight)) %>% ggplot(aes(x=as.numeric(Date), y=Mean, colour=Treatment)) + geom_line() + xlab("Days") + ggtitle("Least Squares Treatment Effect")

pdf("../figs/cattle_example.pdf", width=14)
cattle_spaghetti + cattle_mean
dev.off()

Y <- cattle_data %>% dplyr::select(-c(1, 2)) %>% as.matrix
X <- cattle_data %>% dplyr::select(Treatment) %>% mutate(Treatment = 1*(Treatment == "Treatment")) %>% as.matrix

## Single fit analysis
s <- 2
r <- 9

res <- fit_envelope(Y=Y, X=X, s=s, r=r, maxIters=1000, prior_counts=0, Vinit="OLS", distn="cook")
res$beta_env

tibble(Date=as.numeric(colnames(cattle_data)[-c(1:2)]),
       Envelope=as.numeric(res$beta_env),
       OLS=as.numeric(res$beta_ols)) %>%
    gather(-1, key=Type, value=Coef) %>%
    ggplot(aes(x=Date, y=Coef, col=Type)) + geom_line(size=1) + geom_hline(yintercept=0, linetype="dashed")


envelope_estimator <- tibble(Date=as.numeric(colnames(cattle_data)[-c(1:2)]),
       Control=res$intercept,
       Treatment=as.numeric(res$intercept + res$beta_env)) %>%
    gather(value=Weight, key=Treatment, -Date) %>% ggplot(aes(x=Date, y=Weight, colour=as.factor(Treatment))) + geom_line()

pdf("../figs/envelope_comparison.pdf", width=14)
cattle_mean + envelope_estimator
dev.off()

projected <- as.tibble(Y %*% res$V[, 1:s , drop=FALSE]) %>% bind_cols(as.tibble(X))
projected0 <- as.tibble(Y %*% res$V[, -(1:s) , drop=FALSE]) %>% bind_cols(as.tibble(X))

proj_data <- projected %>% ggplot(aes(x=V1, y=V2, colour=as.factor(Treatment))) + geom_point() + ggtitle("Envelope Subspace")
orth_proj_data <- projected0 %>% ggplot(aes(x=V1, y=V2, colour=as.factor(Treatment))) + geom_point() + ggtitle("Orthogonal Subspace")

pdf("../figs/data_projection.pdf", width=14)
proj_data + orth_proj_data
dev.off()


res2 <- envelopeR::fit_envelope(Y=Y, X=X, s=s, r=r, maxIters=10000, Vinit="OLS", use_py=FALSE)
projected2 <- as.tibble(Y %*% res2$V[, 1:s , drop=FALSE]) %>% bind_cols(as.tibble(X))
proj_data2 <- projected2 %>% ggplot(aes(x=V1, y=V2, colour=as.factor(Treatment))) + geom_point() + ggtitle("Envelope Subspace")
proj_data + proj_data2


projected %>% ggplot(aes(x=V1, y=as.factor(Treatment), fill=as.factor(Treatment))) + geom_density_ridges()

#############################################################################

pdf("../figs/two_week_comparison.pdf")
cattle_data %>% ggplot(aes(x=`84`, y=`98`, colour=Treatment)) + geom_point() +
    ggtitle("Week 8 vs Week 9") + xlab("Week 8") + ylab("Week 9") +
    scale_colour_brewer(type="qual", palette=6)
dev.off()


save_plot("~/Desktop/metab_cartoon.png",
          cattle_data %>% ggplot(aes(x=`84`, y=`98`, colour=Treatment)) +
          geom_point() +
          ggtitle("") +
          xlab("Metabolite 1") + ylab("Metabolite 2") +
          scale_colour_brewer(type="qual", palette=6), base_aspect_ratio=1.3)

pval84 <- t.test(`84` ~ Treatment, cattle_data)$p.value
pval98 <- t.test(`98` ~ Treatment, cattle_data)$p.value
print(pval84)
print(pval98)

pdf("../figs/two_week_mle.pdf", width=14)
cattle_data %>% dplyr::select(Treatment, `84`, `98`) %>% gather(key=Time, value=Weight, -Treatment) %>% ggplot(aes(x=Weight, y=Treatment, fill=Treatment)) + geom_density_ridges() + facet_wrap(~Time)
dev.off()

week8plot <- cattle_data %>% dplyr::select(Treatment, `84`) %>% gather(key=Time, value=Weight, -Treatment) %>% ggplot(aes(x=Weight, y=Treatment, fill=Treatment)) + geom_density_ridges() + theme(legend.pos="None") + ggtitle("Week 8") + ylab("") + xlab("") + scale_fill_brewer(type = "qual", palette = 6, direction = 1)
week9plot <- cattle_data %>% dplyr::select(Treatment, `98`) %>% gather(key=Time, value=Weight, -Treatment) %>% ggplot(aes(x=Weight, y=Treatment, fill=Treatment)) + geom_density_ridges() + theme(legend.pos="None") + ggtitle("Week 9") + ylab("") + xlab("") + scale_fill_brewer(type = "qual", palette = 6, direction = 1)
save_plot("~/Desktop/test.png", plot_grid(week8plot, week9plot, align = "h"), base_aspect_ratio = 2)


Y <- cattle_data %>% dplyr::select(-c(1, 2)) %>% as.matrix
X <- cattle_data %>% dplyr::select(Treatment) %>% mutate(Treatment = 1*(Treatment == "Control")) %>% as.matrix

## Single fit analysis
s <- 1
r <- 1

res <- fit_envelope(Y=Y[, c("84", "98")],
                    X=X,
                    s=s, r=r, maxIters=2000, prior_counts=0)

slp <- res$V[1, 1]/res$V[2, 1]

two_week_envelope_fit <- cattle_data %>% ggplot(aes(x=`84`, y=`98`, colour=Treatment)) + geom_point()  +
    geom_abline(slope=slp, intercept= res$intercept[2] - slp*res$intercept[1], size=1, color="grey", linetype="dashed") + scale_colour_brewer(type="qual", palette=6) +
    xlab("Week 8") + ylab("Week 9") + theme(legend.pos="None")
two_week_envelope_fit


proj_data <- cattle_data %>% dplyr::select(Treatment, `84`, `98`) %>% mutate(Proj = `84` * -res$V[1, 1] - `98` * res$V[1, 2])


two_week_envelope_ridges <- proj_data %>% select(-c(`84`, `98`)) %>%
gather(key=Time, value=Weight, -Treatment) %>% ggplot(aes(x=Weight, y=Treatment, fill=Treatment)) + geom_density_ridges() + scale_fill_brewer(type="qual", palette=6) + theme(legend.pos = "None")

save_plot("~/Desktop/test3.png", plot_grid(two_week_envelope_fit, two_week_envelope_ridges, align="h"), base_aspect_ratio=2)

print(t.test(Proj ~ Treatment, proj_data)$p.value)

cattle_data %>% ggplot(aes(x=`84`, y=`98`, colour=Treatment)) + geom_point()  + geom_abline(slope=res$V[1,1], intercept=res$V[1, 2] + res$intercept[2])

pdf("../figs/two_week_envelope.pdf", width=14)
two_week_envelope_fit + two_week_envelope_ridges
dev.off()

res$beta_env

## bootstrap analysis

fit_resampled <- function(resamp, s, r) {

    df <- resamp$data[resamp$idx, ]

    Y <- df %>% dplyr::select(-(1:2)) %>% as.matrix
    X <- df %>% dplyr::select(Treatment) %>% mutate(Treatment = 1*(Treatment == "Control")) %>% as.matrix
    res <- fit_envelope(Y=Y, X=X, s=s, r=r, maxIters=2000, prior_counts=0)
    res
}

test_out_of_sample <- function(test, params) {

    df <- test$data[test$idx, ]
    Y <- df %>% dplyr::select(-(1:2)) %>% as.matrix
    X <- df %>% dplyr::select(Treatment) %>% mutate(Treatment = 1*(Treatment == "Control")) %>% as.matrix

    beta_ols <- params$beta_ols
    beta_env <- params$beta_env

    Z_ols <- (t(t(Y) - params$intercept) - X %*% beta_ols)
    Z_env <- (t(t(Y) - params$intercept) - X %*% beta_env)

    ## compute MSE
    mse_mle <- mean(Z_ols ^ 2)
    mse_env <- mean(Z_env ^ 2)

    list(
        mse_mle = mse_mle,
        mse_env = mse_env
    )
}

folds <- crossv_kfold(cattle_data, k = 10)

s <- 5
r <- 6
folds %<>% mutate(res = map(train, function(resamp) fit_resampled(resamp, s=s, r=r)))

folds %<>% mutate(prediction = map2(test, res, test_out_of_sample))

folds %<>% mutate(
  mse_env = unlist(map(prediction, function(x) x[[1]])),
  mse_mle = unlist(map(prediction, function(x) x[[2]]))
)

bind_rows(
  summarise(folds, type = "MSE_ENV", mean = mean(mse_env), med = median(mse_env)),
  summarise(folds, type = "MSE_MLE", mean = mean(mse_mle), med = median(mse_mle))
)
