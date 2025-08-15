library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

#' Plot the results of the simulations for the survival resampling inference
#'
#' @param sim_result
#' @param title_mu
#' @param alpha_LCL
#'
#' @returns
#' @export
#'
#' @examples
plot_phase23_results_surv <- function(sim_result, title_mu = expression(lambda[0] == 0 ~ "," ~ lambda[1] == 0 ~ ", "~ lambda[2] == 0), alpha_LCL = 0.025) {
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
  library(cowplot)

  large_text_theme <- theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.1, "cm"),
      legend.key.height = unit(0.6, "cm"),
      plot.title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 16)
    )

  est_by_dose <- sim_result$est_by_dose
  dose_prob_df <- sim_result$dose_prob_df


  # Okabe-Ito color-blind friendly palette (up to 8 colors)
  okabe_ito_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                        "#0072B2", "#D55E00", "#CC79A7", "#999999")

  # Dynamically assign colors to dose levels
  dose_levels <- unique(est_by_dose$dose_selected)
  cb_palette <- setNames(okabe_ito_colors[seq_along(dose_levels)], dose_levels)


  naive_lines <- est_by_dose %>%
    distinct(dose_selected, logHR_naive_bias) %>%
    mutate(linetype = "Naive Estimator")

  bias_plot <- ggplot(est_by_dose, aes(x = w, y = (bias), color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator"), lwd = 1.2) +
    geom_hline(data = naive_lines,
               aes(yintercept = (logHR_naive_bias), color = dose_selected, linetype = linetype), lwd = 1.1) +
    scale_color_manual(name = "Dose Selected", values = cb_palette) +
    scale_linetype_manual(name = "Estimator Type", values = c("Debiased Estimator" = "solid", "Naive Estimator" = "dashed")) +
    ylab("Bias")

  naive_lines_sd <- est_by_dose %>%
    distinct(dose_selected, sd_logHR_naive) %>%
    mutate(linetype = "Naive Estimator")

  min_mse_points <- est_by_dose %>%
    group_by(dose_selected) %>%
    slice_min(mse, with_ties = FALSE)

  naive_lines_mse <- est_by_dose %>%
    distinct(dose_selected, mse_naive) %>%
    mutate(linetype = "Naive Estimator")

  mse_plot <- ggplot(est_by_dose, aes(x = w, y = sqrt(mse), color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator"), lwd = 1.2) +
    geom_hline(data = naive_lines_mse, lwd = 1.1,
               aes(yintercept = sqrt(mse_naive), color = dose_selected, linetype = linetype)) +
    geom_point(data = min_mse_points, aes(x = w, y = sqrt(mse)), shape = 21, fill = "white", size = 3, stroke = 1) +
    geom_text_repel(data = min_mse_points,
                    aes(x = w, y = sqrt(mse), label = paste0("min at w=", round(w, 2))),
                    size = 5,
                    nudge_y = 0.02,
                    direction = "y") +
    scale_color_manual(name = "Dose Selected", values = cb_palette) +
    scale_linetype_manual(values = c("Debiased Estimator" = "solid", "Naive Estimator" = "dashed")) +
    ylab("RMSE") +
    theme_minimal()

  sd_plot <- ggplot(est_by_dose, aes(x = w, y = sd_logHR_tilde, color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator"), lwd = 1.2) +
    geom_hline(data = naive_lines_sd, lwd = 1.1,
               aes(yintercept = sd_logHR_naive, color = dose_selected, linetype = linetype)) +
    scale_color_manual(name = "Dose Selected", values = cb_palette) +
    scale_linetype_manual(values = c("Debiased Estimator" = "solid", "Naive Estimator" = "dashed")) +
    ylab("SD")

  naive_lines_emp_cov <- est_by_dose %>%
    distinct(dose_selected, emp_cov_naive) %>%
    mutate(linetype = "Naive Estimator")

  emp_cov_plot <- ggplot(est_by_dose, aes(x = w, y = emp_cov, color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator")) +
    geom_hline(data = naive_lines_emp_cov,
               aes(yintercept = emp_cov_naive, color = dose_selected, linetype = linetype)) +
    geom_hline(yintercept = 1 - alpha_LCL, linetype = "twodash", color = "blue") +
    scale_color_manual(name = "Dose Selected", values = cb_palette) +
    scale_linetype_manual(values = c("Debiased Estimator" = "solid", "Naive Estimator" = "dashed")) +
    ylab("Empirical Coverage")

  legend <- cowplot::get_legend(bias_plot + large_text_theme + theme(legend.position = "right"))

  bias_plot_clean <- bias_plot + large_text_theme + theme(legend.position = "none")
  sd_plot_clean <- sd_plot + large_text_theme + theme(legend.position = "none")
  mse_plot_clean <- mse_plot + large_text_theme + theme(legend.position = "none")
  emp_cov_plot_clean <- emp_cov_plot + large_text_theme + theme(legend.position = "none")

  plots_combined <- bias_plot_clean / sd_plot_clean / mse_plot_clean +
    plot_annotation(title = title_mu, theme = theme(plot.title = element_text(size = 20, face = "bold")))

  final_plot <- plot_grid(plots_combined, legend, ncol = 2, rel_widths = c(4, 1))

  return(final_plot)
}
