

library(tidyverse)

df <- read_tsv("../data/d_trial_exponential_results.tsv")

# If `asymp_opt_power` is constant per Method, we extract one value per Method
asymp_lines <- df %>%
  group_by(Method) %>%
  summarise(asymp_opt_power = dplyr::first(asymp_opt_power))  # assumes constant per Method

# Main plot
p1 <- ggplot(df, aes(x = d, y = mean, color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.1) +
  # Add horizontal dashed lines for asymptotic power
  geom_hline(data = asymp_lines,
             aes(yintercept = asymp_opt_power, color = Method),
             linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "d",
    y = expression("|" * "<" * hat(v) * "," * v * ">" * "|"^2),
    color = "Method"
  ) +
  theme_bw()


df <- read_tsv("../data/d_trial2_results.tsv")

library(tidyverse)

# If `asymp_opt_power` is constant per Method, we extract one value per Method
asymp_lines <- df %>%
  group_by(Method) %>%
  summarise(asymp_opt_power = dplyr::first(asymp_opt_power))  # assumes constant per Method

# Main plot
p2 <- ggplot(df, aes(x = d, y = mean, color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.1) +
  # Add horizontal dashed lines for asymptotic power
  geom_hline(data = asymp_lines,
             aes(yintercept = asymp_opt_power, color = Method),
             linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "d",
    y = expression("|" * "<" * hat(v) * "," * v * ">" * "|"^2),
    color = "Method"
  ) +
  theme_bw()


library(ggpubr)



#ggarrange(p1, p2, ncol = 2, common.legend =TRUE, legend="bottom")



#Trial 3 performance as a function of d

df <- read_tsv("../data/trial3_results.tsv")

# If `asymp_opt_power` is constant per Method, we extract one value per Method
asymp_lines <- df %>%
  group_by(Method) %>%
  summarise(asymp_opt_power = dplyr::first(asymp_opt_power))  # assumes constant per Method

# Main plot
p3 <- ggplot(df, aes(x = d, y = mean, color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.1) +
  # Add horizontal dashed lines for asymptotic power
  geom_hline(data = asymp_lines,
             aes(yintercept = asymp_opt_power, color = Method),
             linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "d",
    y = expression("|" * "<" * hat(v) * "," * v * ">" * "|"^2),
    color = "Method"
  ) +
  theme_bw()


#p4

#Trial 3 performance as a function of d

df <- read_tsv("../data/trial4_results.tsv")

asymp_lines <- df %>%
  group_by(Method) %>%
  summarise(asymp_opt_power = dplyr::first(asymp_opt_power))

p4 <- ggplot(df, aes(x = M, y = mean, color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.01) +
  geom_line(data = df,
             aes(x = M, y = asymp_opt_power, color = Method),
             linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "M",
    y = expression("|" * "<" * hat(v) * "," * v * ">" * "|"^2),
    color = "Method"
  ) +
  theme_bw()

#Trial 5 performance as a function of M

df <- read_tsv("../data/trial5_results.tsv")

# Main plot
p5 <- ggplot(df, aes(x = M, y = mean, color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.01) +
  # Add horizontal dashed lines for asymptotic power
  geom_line(data = df,
            aes(x = M, y = asymp_opt_power, color = Method),
            linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "M",
    y = expression("|" * "<" * hat(v) * "," * v * ">" * "|"^2),
    color = "Method"
  ) +
  theme_bw()


#Trial 6

df <- read_tsv("../data/trial6_results.tsv")



# If `asymp_opt_power` is constant per Method, we extract one value per Method
asymp_lines <- df %>%
  group_by(Method) %>%
  summarise(asymp_opt_power = dplyr::first(asymp_opt_power))  # assumes constant per Method

# Main plot
p6 <- ggplot(df, aes(x = d, y = mean, color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.1) +
  # Add horizontal dashed lines for asymptotic power
  geom_hline(data = asymp_lines,
             aes(yintercept = asymp_opt_power, color = Method),
             linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "d",
    y = expression("|" * "<" * hat(v) * "," * v * ">" * "|"^2),
    color = "Method"
  ) +
  theme_bw()



#d trial 1 performance as a function of d

df <- read_tsv("../data/d_trial1_results.tsv")


# If `asymp_opt_power` is constant per Method, we extract one value per Method
asymp_lines <- df %>%
  group_by(Method) %>%
  summarise(asymp_opt_power = dplyr::first(asymp_opt_power))  # assumes constant per Method

# Main plot
p7 <- ggplot(df, aes(x = d, y = mean, color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.1) +
  # Add horizontal dashed lines for asymptotic power
  geom_hline(data = asymp_lines,
             aes(yintercept = asymp_opt_power, color = Method),
             linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "d",
    y = expression("|" * "<" * hat(v) * "," * v * ">" * "|"^2),
    color = "Method"
  ) +
  theme_bw()


# Replicate paper plots

ggarrange(p2,p3,
          p6,p5,
          nrow=2,ncol=2,
          common.legend=TRUE,
          legend = "bottom",
          labels=c("a","b","c","d"))



