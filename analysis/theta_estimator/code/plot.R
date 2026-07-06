args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()

data_file <- file.path(script_dir, "..", "data", "theta2_estimation_results.RDS")
combined_plot_file <- file.path(
  script_dir,
  "..",
  "plots",
  "theta2_estimation_combined.pdf"
)

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

results <- readRDS(data_file)

singular_vector_cols <- c(
  "unweighted_stack_svd",
  "estimated_weight_stack_svd",
  "optimal_weight_stack_svd"
)

missing_singular_vector_cols <- setdiff(singular_vector_cols, names(results))
if (length(missing_singular_vector_cols) > 0) {
  stop(
    "Missing singular-vector recovery columns in ",
    data_file,
    ". Rerun theta2_est.R before running plot.R. Missing columns: ",
    paste(missing_singular_vector_cols, collapse = ", ")
  )
}

plot_data <- results |>
  mutate(rel_sq_error = sq_error / theta2^2) |>
  group_by(theta1, theta2) |>
  summarise(
    rrmse = sqrt(mean(rel_sq_error)),
    se = sd(rel_sq_error) / sqrt(n()) / (2 * rrmse),
    .groups = "drop"
  ) |>
  mutate(theta2_label = factor(theta2, levels = c(0.1, 0.5, 0.9)))

p <- ggplot(
  plot_data,
  aes(x = theta1, y = rrmse, color = theta2_label, group = theta2_label)
) +
  geom_errorbar(
    aes(ymin = rrmse - 1.96 * se, ymax = rrmse + 1.96 * se),
    width = 0.015,
    linewidth = 0.7,
    alpha = 0.65
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.2) +
  scale_color_manual(
    name = expression(theta[2]),
    values = c("0.1" = "#0072B2", "0.5" = "#D55E00", "0.9" = "#009E73")
  ) +
  labs(
    x = expression(theta[1]),
    y = expression(RRMSE(hat(theta)[2])),
    title = expression(paste("RRMSE of ", hat(theta)[2], " vs ", theta[1]))
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "gray30"),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p)

singular_vector_data <- results |>
  select(
    theta1,
    theta2,
    all_of(singular_vector_cols)
  ) |>
  pivot_longer(
    cols = all_of(singular_vector_cols),
    names_to = "method",
    values_to = "recovery"
  ) |>
  mutate(
    method = recode(
      method,
      unweighted_stack_svd = "Stack SVD (Unweighted)",
      estimated_weight_stack_svd = "Stack SVD (estimated weights)",
      optimal_weight_stack_svd = "Stack SVD (oracle weights)"
    ),
    method = factor(
      method,
      levels = c(
        "Stack SVD (Unweighted)",
        "Stack SVD (estimated weights)",
        "Stack SVD (oracle weights)"
      )
    ),
    theta2_label = paste0("theta2 = ", theta2)
  ) |>
  group_by(theta1, theta2, theta2_label, method) |>
  summarise(
    mean_recovery = mean(recovery),
    se = sd(recovery) / sqrt(n()),
    .groups = "drop"
  )

p.singular <- ggplot(
  singular_vector_data,
  aes(x = theta1, y = mean_recovery, color = method, group = method)
) +
  geom_errorbar(
    aes(
      ymin = mean_recovery - 1.96 * se,
      ymax = mean_recovery + 1.96 * se
    ),
    width = 0.015,
    linewidth = 0.7,
    alpha = 0.65
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.2) +
  facet_wrap(~ theta2_label, nrow = 1) +
  scale_color_manual(
    name = "Method",
    values = c(
      "Stack SVD (Unweighted)" = "#0072B2",
      "Stack SVD (estimated weights)" = "#D55E00",
      "Stack SVD (oracle weights)" = "#009E73"
    )
  ) +
  labs(
    x = expression(theta[1]),
    y = expression((hat(v)^T * v)^2),
    title = expression(paste("Singular Vector Recovery vs ", theta[1]))
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "gray30"),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 13),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p.singular)

blank_plot <- ggplot() + theme_void()

top_panel <- ggarrange(
  blank_plot,
  p,
  blank_plot,
  nrow = 1,
  ncol = 3,
  widths = c(0.19, 0.62, 0.19)
)

combined_plot <- ggarrange(
  top_panel,
  p.singular,
  nrow = 2,
  ncol = 1,
  labels = c("a", "b"),
  heights = c(1, 1.15)
)

print(combined_plot)

ggsave(combined_plot_file, combined_plot, width = 10, height = 8,
       units = "in", bg = "white")
