setwd(here::here("analysis/r_supp_plots/code"))

M <- 3

c0 <- 1

theta0 <- seq(0.5, 1.5, by=0.0005)

#stacksvd <- ifelse(theta0 > M^{-1/2}*c0^{1/4},
#                   1 - (c0 + M*theta0^2)/(M*theta0^4 + theta0^2),
#                   0)

stacksvd <- ifelse(theta0 > M^{-1/4}*c0^{1/4},
                   1 - (c0 + theta0^2)/(M*theta0^4 + theta0^2),
                   0)

svdstack <- ifelse(theta0 > c0^{1/4},
                   1 - (c0 + theta0^2)/(M*theta0^4 + theta0^2 - (M-1)*c0),
                   0)

df <- data.frame(x=theta0,
                 stacksvd = stacksvd,
                 svdstack = svdstack)

library(reshape2)

df <- melt(df, id.vars = "x")

library(ggplot2)

p <- ggplot(data=df,aes(x=x,y=value, color=variable)) +
  geom_line() +
  theme_bw() +
  geom_vline(xintercept=3^{-1/4}, linetype="dashed", color="grey") +
  geom_vline(xintercept=1, linetype="dashed", color="grey") +
  xlab(expression(theta[0])) +
  ylab("Squared inner product") +
  labs(color="Method")

ggsave(p, filename="../plots/corollary_1_plot.png",
       width=5.71, height=3.85)



