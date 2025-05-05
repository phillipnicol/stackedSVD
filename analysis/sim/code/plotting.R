

df <- read_tsv("../data/trial1_results.tsv")

library(tidyverse)


p <- df |> ggplot(aes(x=M, y = mean, col=Method)) +
  geom_point() + geom_line()






df <- read_tsv("../data/trial3_results.tsv")

library(tidyverse)


p <- df |> ggplot(aes(x=d, y = mean, col=Method)) +
  geom_point() + geom_line() + theme_bw() +
  geom_hline(yintercept = 0.820)






df <- read_tsv("../data/d_trial_exponential_results.tsv")

library(tidyverse)

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

ggarrange(p1, p2, ncol = 2, common.legend =TRUE, legend="bottom")
