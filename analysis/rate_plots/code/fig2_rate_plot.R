
#Panel a

setwd(here::here("analysis/rate_plots/code"))

library(tidyverse)
library(ggpubr)

res <- read.csv("../data/fig2a_rate_overall_results.tsv", sep="\t")

res <- res |> group_by(d, Method) |>
  summarize(mean=mean(Value),
            sd=sd(Value)/sqrt(100),
            asymp_opt_power=mean(asymp_opt_power))


library(tidyverse)

a <- 1.96

p <- ggplot(data=res, aes(x=d, y=mean, color=Method,
                          ymin=mean-a*sd, ymax=mean+a*sd)) +
  geom_point() + geom_line() +
  geom_errorbar(width=0.1) +
  geom_line(aes(x=d, y=asymp_opt_power, color=Method), linetype="dashed") +
  theme_bw() +
  xlab("Dimension (d)") +
  ylab("Squared inner product") +
  scale_x_log10()



# Panel b

res <- read.csv("../data/fig2b_rate_results.tsv", sep="\t")

p2 <- ggplot(data=res, aes(x=d, y=mean, color=Method,
                           ymin = mean-a*std_err,
                           ymax=mean+a*std_err)) +
  geom_point() + geom_line() +
  geom_errorbar(width=0.1) +
  geom_line(aes(x=d, y=asymp_opt_power, color=Method), linetype="dashed") +
  theme_bw() +
  xlab("Dimension (d)") +
  ylab("Squared inner product") +
  scale_x_log10()

res <- read.csv("../data/fig2c_rate_results.tsv", sep="\t")




p3 <- ggplot(data=res, aes(x=d, y=mean, color=Method,
                           ymin = mean-a*std_err,
                           ymax=mean+a*std_err)) +
  geom_point() + geom_line() +
  geom_errorbar(width=0.1) +
  geom_line(aes(x=d, y=asymp_opt_power, color=Method), linetype="dashed") +
  theme_bw() +
  xlab("Dimension (d)") +
  ylab("Squared inner product") +
  scale_x_log10()



library(ggpubr)

p.full <- ggarrange(p, p2, p3, nrow=1, ncol=3, common.legend=TRUE,
                    labels=c("a", "b", "c"))


ggsave(p.full, filename="../plots/fig2_rversion.png",
       width=11.3, height=3.95)

ggsave(p.full, filename="../plots/fig2_rversion.pdf",
       width=11.3, height=3.95)
