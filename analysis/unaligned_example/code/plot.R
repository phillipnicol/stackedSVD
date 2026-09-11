args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()

data_dir <- file.path(script_dir, "..", "data")
plot_dir <- file.path(script_dir, "..", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)

color_map <- c(
  "Weighted Stack-SVD" = "#aec7e8",
  "Unweighted Stack-SVD" = "#1f77b4",
  "Weighted SVD-Stack" = "#ffbb78",
  "Unweighted SVD-Stack" = "#ff7f0e"
)

n_empirical_points <- 16
base_c <- 1
kink_colors <- c(
  s_opt = "#756bb1",
  s_dagger = "#31a354"
)

read_theta_results <- function(theta_values) {
  result_files <- file.path(
    data_dir,
    paste0("unaligned_simulation_results_theta_", theta_values, ".csv")
  )

  missing_files <- result_files[!file.exists(result_files)]
  if (length(missing_files) > 0) {
    stop("Missing expected theta result CSV(s): ", paste(missing_files, collapse = ", "))
  }

  lapply(seq_along(result_files), function(i) {
    read.csv(result_files[i]) |>
      mutate(
        theta = theta_values[i],
        theta_label = paste0("theta = ", theta_values[i])
      )
  }) |>
    bind_rows()
}

prepare_plot_data <- function(results) {
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

  empirical_psi_grid <- summary_stats |>
    distinct(theta, psi, sin_psi) |>
    group_by(theta) |>
    arrange(sin_psi, .by_group = TRUE) |>
    mutate(
      psi_index = row_number(),
      keep_point = psi_index %in% unique(round(seq(
        1,
        n(),
        length.out = min(n_empirical_points, n())
      )))
    ) |>
    ungroup() |>
    filter(keep_point) |>
    select(theta, psi)

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
    inner_join(empirical_psi_grid, by = c("theta", "psi")) |>
    select(theta, theta_label, psi, sin_psi, method, mean_result, sd_result) |>
    mutate(method = ifelse(
      method == "stackedSVD_mean",
      "Unweighted Stack-SVD",
      "Unweighted SVD-Stack"
    ))

  theoretical_data <- summary_stats |>
    pivot_longer(
      cols = c(stackedSVD_theoretical, SVDstacked_theoretical),
      names_to = "method",
      values_to = "theoretical_value"
    ) |>
    select(theta, theta_label, psi, sin_psi, method, theoretical_value) |>
    mutate(method = ifelse(
      method == "stackedSVD_theoretical",
      "Unweighted Stack-SVD",
      "Unweighted SVD-Stack"
    ))

  list(empirical = empirical_data, theoretical = theoretical_data)
}

stacksvd_kinks <- function(theta_values) {
  data.frame(theta = theta_values) |>
    mutate(
      theta_label = paste0("theta = ", theta),
      theta4 = theta^4,
      kink_type = case_when(
        theta4 > 2 * base_c ~ "s_opt",
        theta4 > base_c / 2 & theta4 < 2 * base_c ~ "s_dagger",
        TRUE ~ NA_character_
      ),
      kink_s = case_when(
        theta4 > 2 * base_c ~ 1 - sqrt(2 * base_c) / theta^2,
        theta4 > base_c / 2 & theta4 < 2 * base_c ~
          sqrt(2 * base_c) / theta^2 - 1,
        TRUE ~ NA_real_
      ),
      kink_label = case_when(
        kink_type == "s_opt" ~ "s['*']",
        kink_type == "s_dagger" ~ "s['†']",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(kink_s), kink_s >= 0, kink_s <= 1)
}

make_unaligned_plot <- function(theta_values, facet = TRUE) {
  plot_data <- read_theta_results(theta_values) |>
    prepare_plot_data()
  kink_data <- stacksvd_kinks(theta_values)

  p <- ggplot() +
    geom_errorbar(
      data = plot_data$empirical,
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
      data = plot_data$empirical,
      aes(x = sin_psi, y = mean_result, color = method),
      size = 3,
      alpha = 0.9
    ) +
    geom_line(
      data = plot_data$theoretical,
      aes(x = sin_psi, y = theoretical_value, linetype = method),
      linewidth = 1,
      color = "black",
      alpha = 0.6
    ) +
    geom_vline(
      data = filter(kink_data, kink_type == "s_opt"),
      aes(xintercept = kink_s),
      linetype = "dashed",
      linewidth = 0.5,
      color = kink_colors[["s_opt"]],
      alpha = 0.35
    ) +
    geom_text(
      data = filter(kink_data, kink_type == "s_opt"),
      aes(x = kink_s, y = Inf, label = kink_label),
      parse = TRUE,
      color = kink_colors[["s_opt"]],
      angle = 90,
      hjust = 1.15,
      vjust = -0.35,
      size = 3.4,
      alpha = 0.8
    ) +
    geom_vline(
      data = filter(kink_data, kink_type == "s_dagger"),
      aes(xintercept = kink_s),
      linetype = "dashed",
      linewidth = 0.5,
      color = kink_colors[["s_dagger"]],
      alpha = 0.35
    ) +
    geom_text(
      data = filter(kink_data, kink_type == "s_dagger"),
      aes(x = kink_s, y = Inf, label = kink_label),
      parse = TRUE,
      color = kink_colors[["s_dagger"]],
      angle = 90,
      hjust = 1.15,
      vjust = -0.35,
      size = 3.4,
      alpha = 0.8
    ) +
    scale_linetype_manual(
      name = "Theoretical",
      values = c(
        "Unweighted Stack-SVD" = "dashed",
        "Unweighted SVD-Stack" = "dashed"
      ),
      guide = "none"
    ) +
    scale_color_manual(
      name = "Method",
      values = color_map[c("Unweighted Stack-SVD", "Unweighted SVD-Stack")]
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

  if (facet) {
    p <- p + facet_wrap(~ theta_label, nrow = 1, scales = "free_y")
  } else {
    p <- p + ggtitle(unique(plot_data$empirical$theta_label)) +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  }

  p
}

theta_single <- 0.9
theta_panels <- c(1.15, 1.5, 3)

p_single <- make_unaligned_plot(theta_single, facet = FALSE)
print(p_single)

ggsave(
  file.path(plot_dir, "unaligned_simulation_comparison_theta_0.9.pdf"),
  p_single,
  width = 5.2,
  height = 3.95,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "unaligned_simulation_comparison_theta_0.9.png"),
  p_single,
  width = 5.2,
  height = 3.95,
  dpi = 300,
  bg = "white"
)

p_panels <- make_unaligned_plot(theta_panels, facet = TRUE)
print(p_panels)

ggsave(
  file.path(plot_dir, "unaligned_simulation_comparison_theta_panels.pdf"),
  p_panels,
  width = 11.3,
  height = 3.95,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "unaligned_simulation_comparison_theta_panels.png"),
  p_panels,
  width = 11.3,
  height = 3.95,
  dpi = 300,
  bg = "white"
)
