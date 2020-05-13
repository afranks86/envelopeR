library(tidyverse)

load("../results/scores100_complete.Rdata")

score_tib  <- as_tibble(reshape2::melt(score_array))

colnames(score_tib)[1:3]  <- c("NumCov", "Rep", "Beta")


boot_sd <- function(val, nboot=1000) {
  sd(sapply(1:nboot, function(i) mean(sample(val, length(val), replace=TRUE), na.rm=TRUE)))
}


score_tib  %>% group_by(NumCov, Beta)  %>%
  summarize(mean_val = mean(value, na.rm=TRUE), sd_val=boot_sd(value), max_val = max(value, na.rm=TRUE))  %>%
  ungroup() %>%
  ggplot(aes(x=Beta, col=as.factor(NumCov))) + geom_line(aes(y=mean_val)) +
  geom_errorbar(aes(ymin=mean_val-1.96*sd_val,
                    ymax=mean_val+1.96*sd_val),
                width=.2) + theme_bw(base_size=16) +
  ylim(c(0, 1)) +
  ylab("Subspace Similarity") +
  xlab("Signal Strength of the Mean Coefficients") +
  scale_x_continuous() +
  guides(col=guide_legend(title="Number of\nCovariates")) ->
  scores_plot
scores_plot
ggsave(scores_plot, file="../figs/scores_plot.pdf")


##  geom_point(aes(y=max_val, col=as.factor(NumCov))) +


## score_array[2, 70:100, ] <- score_mat2


load("../results/remote/sim_results2.Rdata")

make_plot  <- function(type, array, dnames=NULL) {

  if(is.null(dnames)) {
    dnames  <- list(Method = c("Ignoring Mean", "Envelope"),
                    NumCov = 1:4,
                    Rep = 1:100,
                    BetaSd = 1:5)
  }
  if(type %in% c("stein", "squared_error")) {

    ## if(type == "stein")
    ##   array  <- steins_loss_array
    ## else if(type == "squared_error")
    ##   array  <- squared_error_loss_array

    dimnames(array)  <- dnames

    array_long  <- as_tibble(as.tbl_cube(array))
    array_long  <- array_long %>% rename(loss = array )

    loss_plot <- array_long %>%
      pivot_wider(names_from = "Method", values_from = "loss") %>%
      mutate(Ratio = `Ignoring Mean`/Envelope) %>%
      group_by(NumCov, BetaSd) %>%
      summarize(mean_loss = mean(Ratio, na.rm=TRUE),
                med_loss = median(Ratio, na.rm=TRUE)) %>%
      ggplot() + geom_line(aes(x=BetaSd, y=mean_loss, col=as.factor(NumCov))) +
      theme_bw()
    
  } else {
    array  <- subspace_sim_array

    dimnames(array)  <- dnames

    array_long  <- as_tibble(as.tbl_cube(array))
    array_long  <- array_long %>% rename(loss = array )

    array_long  %>% group_by(Method, NumCov, BetaSd)  %>%
      summarize(mean_val = mean(loss, na.rm=TRUE), sd_val=boot_sd(loss), max_val = max(loss, na.rm=TRUE))  %>%
      ungroup() %>%
      ggplot(aes(x=BetaSd, col=as.factor(NumCov))) + geom_line(aes(y=mean_val, lty=Method)) +
      ## geom_errorbar(aes(ymin=mean_val-1.96*sd_val,
      ##                   ymax=mean_val+1.96*sd_val),
      ##               width=.2) +
      theme_bw(base_size=16) +
      ylim(c(0, 1)) +
      ylab("Subspace Similarity") +
      xlab("Signal Strength of the Mean Coefficients") +
      scale_x_continuous() +
      guides(col=guide_legend(title="Number of\nCovariates")) ->
      loss_plot



  }

  loss_plot
}

make_plot("stein")
make_plot("squared_error")
make_plot("subspace")

subspace_sim_array[, 1, 1:44, 5]
subspace_sim_array[, 2, 1:44, 5]


### Comparison to full join model
load("../results/remote/sim_s_comparison.Rdata")
### SECTION 2 SIM RESULTS

array  <- steins_loss_array[, , , 1]
dnames  <- list(Stilde = c(4, 8, 12, 25),
                Rep = 1:100,
                NumCov = 1:4)
dimnames(array) <-  dnames

array_long  <- as_tibble(as.tbl_cube(array))
array_long  <- array_long %>% rename(loss = array )
loss_plot <- array_long %>%
  pivot_wider(names_from = "Stilde", values_from = "loss") %>%
  mutate_at(-(1:2), function(x) x - .$`4`) %>%
  group_by(NumCov) %>%
  summarize_at(-1, list(function(x) mean(x, na.rm=TRUE), boot_sd)) %>%
  pivot_longer(cols=-1,
               names_to = c("Stilde", ".value"),
               names_pattern = "([0-9]+)_fn(.)",
               values_to="Ratio") %>%
  ggplot(aes(x=as.numeric(Stilde), y=`1`)) +
  geom_line(aes(col=as.factor(NumCov)), size=1.5) +
  geom_errorbar(aes(ymin=`1`-1.96*`2`,
                    ymax=`1`+1.96*`2`),
                width=.2) + theme_bw(base_size=16) +
  labs(colour="q") +
  ylab("Percent increase in risk") +
  xlab(expression(tilde(s)))
loss_plot
ggsave(file="../figs/s_misspecification.pdf", loss_plot)



#################################


dim(squared_error_loss_array)

apply(squared_error_loss_array[, , 1, ], 1, function(x) mean(x, na.rm=TRUE))
apply(squared_error_loss_array[, , 2, ], 1, function(x) mean(x, na.rm=TRUE))
apply(squared_error_loss_array[, , 3, ], 1, function(x) mean(x, na.rm=TRUE))

apply(steins_loss_array[, , 1, ], 1, function(x) median(x, na.rm=TRUE))


steins1  <- t(steins_loss_array[, , 2, ])
colnames(steins1)  <- c(4, 8, 12, 25)


steins1 %>% as_tibble %>% mutate_all(., function(x) (x - .$`4`)/.$`4`) %>%
  pivot_longer(everything()) %>%
  mutate(name = fct_reorder(name, as.numeric(name))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_violin() + theme_bw()

se1  <- t(squared_error_loss_array[, , 1, ])
colnames(se1)  <- c(4, 8, 12, 25)
se1 %>% as_tibble %>% mutate_all(., function(x) (x - .$`4`)/.$`4`) %>%
  pivot_longer(everything()) %>%
  mutate(name = fct_reorder(name, as.numeric(name))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_violin() + theme_bw()




library(ggridges)
as_tibble(steins1) %>% mutate_all(., function(x) x - .$`4`) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(y=name, x=value)) +
  geom_density_ridges() +
  theme_bw()


#################################################
load("../results/remote/sim_mean_results.Rdata")

array  <- steins_loss_array[, , , 1]
dnames  <- list(Method = c("Ignoring Mean", "Envelope"),
                    NumCov = 1:4,
                    Rep = 1:100)

dimnames(array) <-  dnames
array_long  <- as_tibble(as.tbl_cube(array)) %>% rename(loss = array)

array_long %>% pivot_wider(names_from = "Method", values_from = loss) %>%
  mutate(Percent = (`Ignoring Mean` - Envelope)/Envelope) %>%
  ggplot() + geom_boxplot(aes(x=as.factor(NumCov), y=Percent)) + theme_bw() +
  ylim(c(0, 2)) ->
loss_plot

## loss_plot <- array_long %>%
##   pivot_wider(names_from = "Method", values_from = "loss") %>%
##   mutate(Ratio = `Ignoring Mean`/Envelope) %>%
##   group_by(BetaSd) %>%
##   summarize(mean_loss = mean(Ratio, na.rm=TRUE),
##             med_loss = median(Ratio, na.rm=TRUE)) %>%
##   ggplot() + geom_line(aes(x=BetaSd, y=mean_loss)) +
##   theme_bw()

## loss_plot

####################################


apply(steins_loss_array[, 1, , ], c(1, 3), function(x) mean(x, na.rm=TRUE))

mean((steins_loss_array[1, 1, , 1] - steins_loss_array[2, 1, , 1]) /steins_loss_array[2, 1, , 1], na.rm=TRUE)


array  <- steins_loss_array[, , , 1]


array_long  <- as_tibble(as.tbl_cube(array))
array_long  <- array_long %>% rename(loss = array )



make_plot("steins", array, dnames)

dimnames(array)  <- dnames



array_long  %>%
  pivot_wider(names_from=Method, values_from=loss) %>%
  mutate(Ratio = `Ignoring Mean`/Envelope) %>%
  group_by(BetaSd)  %>%
  summarize(mean_val = mean(Ratio, na.rm=TRUE), sd_val=boot_sd(Ratio),
            max_val = max(Ratio, na.rm=TRUE))  %>%
  ungroup() %>%
  ggplot(aes(x=BetaSd)) + geom_line(aes(y=mean_val)) +
  geom_errorbar(aes(ymin=mean_val-1.96*sd_val,
                    ymax=mean_val+1.96*sd_val),
                width=.2) + theme_bw(base_size=16) +
  ylab("Subspace Similarity") +
  xlab("Signal Strength of the Mean Coefficients") +
  scale_x_continuous() +
  guides(col=guide_legend(title="Number of\nCovariates")) ->
  scores_plot
