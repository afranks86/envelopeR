library(tidyverse)
library(cubelyr)

boot_sd <- function(val, nboot=1000) {
  sd(sapply(1:nboot, function(i) mean(sample(val, length(val), replace=TRUE), na.rm=TRUE)))
}



### Comparison to full join model
## load("../results/remote/sim_s_comparison.Rdata")
load("../results/sim_s_comparison_2020-10-16.Rdata")
### SECTION 2 SIM RESULTS

array  <- steins_loss_array[, , , 1]
dnames  <- list(Stilde = c(2, 3, 4, 8, 12, 25),
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
ggsave(file="../figs/s_misspecification_full.pdf", loss_plot)


## Goodnness of fit test
load("../results/sim_s_comparison_2020-10-16.Rdata")



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

####################################

# Sensitivity to initial values

library(tidyverse)
library(cubelyr)
load("../results/init_sims_full.Rdata")

# pars <- sapply(params, function(x) bquote(beta ~ .(x["beta"]) ~ gamma == .(x["gamma"])))
pars <- sapply(params, function(x) sprintf("\u03b2=%.2f, \u03a9=%.2f", round(x["beta"], 2), round(x["gamma"], 2)))

dimnames(subspace_init_mat) <- dimnames(subspace_sim_mat) <- list(Param=pars, Type=c("Random", "OLS"), Rep=paste0("Rep", 1:10))

sim_tib <- as.tbl_cube(subspace_sim_mat, met="Score") %>% as_tibble
sim_tib$Class <- "Final"
init_tib <- as.tbl_cube(subspace_init_mat, met="Score") %>% as_tibble
init_tib$Class <- "Initial"
rbind(sim_tib, init_tib) %>%
  ggplot(aes(x=fct_relevel(Class, "Initial", "Final"), y=Score, color=Type, group=interaction(Rep, Type))) +
  geom_line(position=position_jitter(w=0.1, h=0)) + theme_minimal(base_size=16) +
  ## geom_point(size=0.5) +
  facet_wrap(~ Param) +
  xlab("") + ggtitle("Importance of Initialization") ->
  init_plot

ggsave(init_plot, file="../figs/initialization.pdf")

####################################
# Run time results

library(tidyverse)
library(cubelyr)
load("../results/benchmark_times.RData")

benchmark_tib <- as.tbl_cube(benchmark_array, met="time") %>% as_tibble

benchmark_tib %>% ggplot(aes(x=p, y=time, col=as.factor(s))) +
  geom_point(position=position_jitter(h=0)) +
  theme_bw(base_size=16) +
  ylab("time (seconds)") +
  labs(col="s") ->
  time_benchmark


ggsave(init_plot, file="../figs/initialization.pdf")
