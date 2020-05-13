library(tidyverse)




boot_sd <- function(val, nboot=1000) {
  sd(sapply(1:nboot, function(i) mean(sample(val, length(val), replace=TRUE), na.rm=TRUE)))
}



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


#################################################

load("../results/remote/sim_mean_results.Rdata")

array  <- steins_loss_array[, , , 5]
dnames  <- list(Method = c("Ignoring Mean", "Envelope"),
                    NumCov = 1:4,
                    Rep = 1:100)

dimnames(array) <-  dnames
array_long  <- as_tibble(as.tbl_cube(array)) %>% rename(loss = array)

array_long %>% pivot_wider(names_from = "Method", values_from = loss) %>%
  mutate(Percent = (`Ignoring Mean` - Envelope)/Envelope) %>%
  ggplot() + geom_boxplot(aes(x=as.factor(NumCov), y=Percent, fill=as.factor(NumCov))) + theme_bw(base_size=16) +
  xlab("q") + ylab("Percent Increase in Risk") +
  labs(fill="q") +
  #ylim(c(0, 2)) ->
  loss_plot

loss_plot
ggsave(file="../figs/mean_results.pdf", loss_plot, width=6)

####################################


load("../results/scores100_complete.Rdata")

score_tib  <- as_tibble(reshape2::melt(score_array))

colnames(score_tib)[1:3]  <- c("NumCov", "Rep", "Beta")


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

