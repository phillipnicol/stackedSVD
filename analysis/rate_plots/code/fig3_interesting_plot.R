
#Panel a

res <- read.csv("../data/fig3a_interesting_results.tsv", sep="\t")

library(tidyverse)

a <- 1.96

p <- ggplot(data=res, aes(x=M, y=mean, color=Method,
                          ymin = mean-a*std_err,
                          ymax=mean+a*std_err)) +
  geom_point() + geom_line() +
  geom_errorbar(width=0.66) +
  geom_line(aes(x=M, y=asymp_opt_power, color=Method), linetype="dashed") +
  theme_bw() +
  xlab("Number of matrices (M)") +
  ylab("Squared inner product")



# Panel b

res <- read.csv("../data/fig3b_interesting_results.tsv", sep="\t")

p2 <- ggplot(data=res, aes(x=M, y=mean, color=Method,
                           ymin = mean-a*std_err,
                           ymax=mean+a*std_err)) +
  geom_point() + geom_line() +
  geom_errorbar(width=0.66) +
  geom_line(aes(x=M, y=asymp_opt_power, color=Method), linetype="dashed") +
  theme_bw() +
  xlab("Number of matrices (M)") +
  ylab("Squared inner product")

res <- read.csv("../data/fig3c_interesting_results.tsv", sep="\t")




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


ggsave(p.full, filename="../plots/fig3_rversion.png",
       width=11.3, height=3.95)

ggsave(p.full, filename="../plots/fig3_rversion.pdf",
       width=11.3, height=3.95)
