setwd(here::here("analysis/rate_plots/code"))

library(tidyverse)
library(ggpubr)

a <- 1.96

color_map <- c(
  "Weighted Stack-SVD" = "#aec7e8",
  "Unweighted Stack-SVD" = "#1f77b4",
  "Weighted SVD-Stack" = "#ffbb78",
  "Unweighted SVD-Stack" = "#ff7f0e"
)

# Panel a: fixed d, increasing M
res <- read.csv("../data/trial2_exp_results.tsv", sep="\t")

p <- ggplot(data=res, aes(x=M, y=mean, color=Method,
                          ymin=mean-a*std_err,
                          ymax=mean+a*std_err)) +
  geom_point() +
  geom_line() +
  geom_errorbar(width=0.66) +
  geom_line(aes(x=M, y=asymp_opt_power, color=Method), linetype="dashed") +
  theme_bw() +
  scale_color_manual(values=color_map) +
  xlab("Number of matrices (M)") +
  ylab("Squared inner product")

# Panel b: fixed M, increasing d
res <- read.csv("../data/d_trial_exponential_results.tsv", sep="\t")

p2 <- ggplot(data=res, aes(x=d, y=mean, color=Method,
                           ymin=mean-a*std_err,
                           ymax=mean+a*std_err)) +
  geom_point() +
  geom_line() +
  geom_errorbar(width=0.1) +
  geom_line(aes(x=d, y=asymp_opt_power, color=Method), linetype="dashed") +
  theme_bw() +
  scale_color_manual(values=color_map) +
  xlab("Dimension (d)") +
  ylab("Squared inner product") +
  scale_x_log10()

p.full <- ggarrange(p, p2, nrow=1, ncol=2, common.legend=TRUE,
                    labels=c("a", "b"))

ggsave(p.full, filename="../plots/exponential_rate_plot.png",
       width=8.0, height=3.95)

ggsave(p.full, filename="../plots/exponential_rate_plot.pdf",
       width=8.0, height=3.95)
