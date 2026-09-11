setwd(here::here("analysis/sim/code"))

library(tidyverse)
library(ggpubr)
library(magick)
library(grid)

method_colors <- c(
  "Stack SVD" = "#1f77b4",
  "SVD Stack" = "#ff7f0e",
  "Weighted Stack SVD" = "#aec7e8",
  "Weighted SVD Stack" = "#ffbb78"
)

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

df$Var2 <- factor(as.character(df$Var2), levels=names(method_colors))

df <- df |> filter(sd == 0.1) #Try filtering to sd = 1

df$mean <- paste("Mean = ", round(log(df$mean), digits=1))
df$sd <- paste("Sd = ", df$sd)

df$mean_lab <- factor(df$mean, labels = paste0("mu == ", levels(df$mean)))
df$sd_lab <- factor(df$sd, labels = paste0("sigma == ", levels(df$sd)))


p <- df |> ggplot(aes(x=Var2, y=value, fill=Var2)) +
  #facet_grid(mean~sd,
  #           scales="free_y") +
  facet_wrap(~mean) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color="grey50") +
  theme_bw() +
  scale_fill_manual(values=method_colors) +
  ylab("Bias") +
  xlab("") +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





#ggsave(p, filename="../plots/asymptotic_bias.png",
#       width=8.46, height=7.44, units="in")













df <- readRDS("../data/ctheta_sim_var_1.1_0.1.RDS")
df$mean <- 1
df$sd <- 0.1

df1 <- readRDS("../data/ctheta_sim_var_1.1_0.25.RDS")
df1$mean <- 1
df1$sd <- 0.25

df2 <- readRDS("../data/ctheta_sim_var_1.1_1.RDS")
df2$mean <- 1
df2$sd <- 1

df3 <- readRDS("../data/ctheta_sim_var_1.5_0.1.RDS")
df3$mean <- 1.5
df3$sd <- 0.1


df4 <- readRDS("../data/ctheta_sim_var_1.5_0.25.RDS")
df4$mean <- 1.5
df4$sd <- 0.25

#df5 <- readRDS("../data/ctheta_sim_var_1.5_1.RDS")
#df5$mean <- 1.5
#df5$sd <- 1

df6 <- readRDS("../data/ctheta_sim_var_2_0.1.RDS")
df6$mean <- 2
df6$sd <- 0.1

df7 <- readRDS("../data/ctheta_sim_var_2_0.25.RDS")
df7$mean <- 2
df7$sd <- 0.25

#df8 <- readRDS("../data/ctheta_sim_var_2_1.RDS")
#df8$mean <- 2
#df8$sd <- 1


#Bind all of the data frames into one df
#df <- rbind(df, df1, df2, df3, df4, df5, df6, df7, df8)
df <- rbind(df, df1, df2, df3, df4, df6, df7)

df$Var2 <- factor(c("Stack SVD", "SVD Stack", "Weighted Stack SVD", "Weighted SVD Stack")[df$Var2],
                  levels=names(method_colors))


df <- df |> filter(sd == 0.1) #Try filtering to sd = 1

df$mean <- paste("Mean = ", round(log(df$mean), digits=1))
df$sd <- paste("Sd = ", df$sd)

df$mean_lab <- factor(df$mean, labels = paste0("mu == ", levels(df$mean)))
df$sd_lab <- factor(df$sd, labels = paste0("sigma == ", levels(df$sd)))


p.var <- df |> ggplot(aes(x=Var2, y=sqrt(value), fill=Var2)) +
  #facet_grid(mean~sd,
  #           scales="free_y") +
  facet_wrap(~mean) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color="grey50") +
  theme_bw() +
  scale_fill_manual(values=method_colors) +
  ylab("Standard error") +
  xlab("") +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




legend_img <- image_read_pdf("../plots/fig3_legend_v3.pdf", density=300) |>
  image_background("white", flatten=TRUE) |>
  image_trim(fuzz=2)
legend_grob <- rasterGrob(as.raster(legend_img), interpolate=TRUE)
legend_plot <- as_ggplot(legend_grob)

p.panels <- ggarrange(p, p.var, nrow=2, labels=c("a","b"))

p.big <- ggarrange(p.panels, legend_plot, nrow=1, widths=c(1, 0.28))


ggsave(p.big, filename="../plots/asymptotic_bias.png",
       width=9.6, height=7.34, units="in", bg="white")


ggsave(p.big, filename="../plots/asymptotic_bias.pdf",
       width=9.6, height=7.34, units="in", bg="white")
