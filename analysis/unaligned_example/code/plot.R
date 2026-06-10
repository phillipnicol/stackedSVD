args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()

data_dir <- file.path(script_dir, "..", "data")
plot_dir <- file.path(script_dir, "..", "plots")

library(ggplot2)
library(dplyr)
library(tidyr)

theta_values <- c(1.25, 1.75, 3)
result_files <- file.path(
  data_dir,
  paste0("unaligned_simulation_results_theta_", theta_values, ".csv")
)

missing_files <- result_files[!file.exists(result_files)]
if (length(missing_files) > 0) {
  stop("Missing expected theta result CSV(s): ", paste(missing_files, collapse = ", "))
}

results <- lapply(seq_along(result_files), function(i) {
  read.csv(result_files[i]) |>
    mutate(
      theta = theta_values[i],
      theta_label = paste0(theta_values[i], "")
    )
}) |>
  bind_rows()

summary_stats <- results |>
  group_by(theta, theta_label, psi) |>
  summarise(
    stackedSVD_mean = mean(stackedSVD_result),
    stackedSVD_sd = sd(stackedSVD_result),
    SVDstacked_mean = mean(SVDstacked_result),
    SVDstacked_sd = sd(SVDstacked_result),
    stackedSVD_theoretical = first(stackedSVD_theoretical),
    SVDstacked_theoretical = first(SVDstacked_theoretical),
    .groups = "drop"
  ) |>
  mutate(sin_psi = sin(psi))

empirical_data <- summary_stats |>
  pivot_longer(
    cols = c(stackedSVD_mean, SVDstacked_mean),
    names_to = "method",
    values_to = "mean_result"
  ) |>
  pivot_longer(
    cols = c(stackedSVD_sd, SVDstacked_sd),
    names_to = "method_sd",
    values_to = "sd_result"
  ) |>
  filter(
    (method == "stackedSVD_mean" & method_sd == "stackedSVD_sd") |
      (method == "SVDstacked_mean" & method_sd == "SVDstacked_sd")
  ) |>
  select(theta, theta_label, psi, sin_psi, method, mean_result, sd_result) |>
  mutate(
    method = ifelse(method == "stackedSVD_mean", "Stack SVD", "SVD Stack"),
    theta_label = paste0("theta = ", theta_label)
  )

theoretical_data <- summary_stats |>
  pivot_longer(
    cols = c(stackedSVD_theoretical, SVDstacked_theoretical),
    names_to = "method",
    values_to = "theoretical_value"
  ) |>
  select(theta, theta_label, psi, sin_psi, method, theoretical_value) |>
  mutate(
    method = ifelse(method == "stackedSVD_theoretical", "Stack SVD", "SVD Stack"),
    theta_label = paste0("theta = ", theta_label)
  )

p <- ggplot() +
  geom_errorbar(
    data = empirical_data,
    aes(
      x = sin_psi,
      ymin = mean_result - sd_result,
      ymax = mean_result + sd_result,
      color = method
    ),
    width = 0.025,
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_point(
    data = empirical_data,
    aes(x = sin_psi, y = mean_result, color = method),
    size = 3,
    alpha = 0.9
  ) +
  geom_line(
    data = theoretical_data,
    aes(x = sin_psi, y = theoretical_value, linetype = method),
    linewidth = 1,
    color = "black",
    alpha = 0.6
  ) +
  facet_wrap(~ theta_label, nrow = 1, scales = "free_y") +
  scale_linetype_manual(
    name = "Theoretical",
    values = c("Stack SVD" = "dashed", "SVD Stack" = "dashed"),
    guide = "none"
  ) +
  scale_color_manual(
    name = "Method",
    values = c("Stack SVD" = "#1b9e77", "SVD Stack" = "#d95f02")
  ) +
  labs(
    x = expression(sin(psi)),
    y = "Squared Frobenius Norm"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold")
  )

print(p)

ggsave(file.path(plot_dir, "unaligned_simulation_comparison_theta_panels.pdf"),
       p, width = 11.3, height = 3.95, dpi = 300, bg = "white")

ggsave(file.path(plot_dir, "unaligned_simulation_comparison_theta_panels.png"),
       p, width = 11.3, height = 3.95, dpi = 300, bg = "white")
