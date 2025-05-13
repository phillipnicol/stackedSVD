df <- readRDS("../data/ctheta_sim_1.1_0.1.RDS")
df$mean <- 1
df$sd <- 0.1

df1 <- readRDS("../data/ctheta_sim_1.1_0.25.RDS")
df1$mean <- 1
df1$sd <- 0.25

df2 <- readRDS("../data/ctheta_sim_1.1_1.RDS")
df2$mean <- 1
df2$sd <- 1

df3 <- readRDS("../data/ctheta_sim_1.5_0.1.RDS")
df3$mean <- 1.5
df3$sd <- 0.1


df4 <- readRDS("../data/ctheta_sim_1.5_0.25.RDS")
df4$mean <- 1.5
df4$sd <- 0.25

df5 <- readRDS("../data/ctheta_sim_1.5_1.RDS")
df5$mean <- 1.5
df5$sd <- 1

df6 <- readRDS("../data/ctheta_sim_2_0.1.RDS")
df6$mean <- 2
df6$sd <- 0.1

df7 <- readRDS("../data/ctheta_sim_2_0.25.RDS")
df7$mean <- 2
df7$sd <- 0.25

df8 <- readRDS("../data/ctheta_sim_2_1.RDS")
df8$mean <- 2
df8$sd <- 1


#Bind all of the data frames into one df
df <- rbind(df, df1, df2, df3, df4, df5, df6, df7, df8)

df$mean <- paste("Mean = ", df$mean)
df$sd <- paste("Sd = ", df$sd)

df$mean_lab <- factor(df$mean, labels = paste0("mu == ", levels(df$mean)))
df$sd_lab <- factor(df$sd, labels = paste0("sigma == ", levels(df$sd)))



library(tidyverse)

p <- df |> ggplot(aes(x=Var2, y=value)) +
  facet_grid(mean~sd) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color="blue") +
  theme_bw() +
  ylab("Bias") +
  xlab("")

#Rotate x-axis labels 45 deg
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
