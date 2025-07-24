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

  est_by_dose <- sim_result$est_by_dose
  dose_prob_df <- sim_result$dose_prob_df

  # Prepare the data for the horizontal lines
  naive_lines <- est_by_dose %>%
    distinct(dose_selected, logHR_naive_bias) %>%
    mutate(linetype = "Naive Estimator")

  # Bias plot
  bias_plot <- ggplot(est_by_dose, aes(x = w, y = bias, color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator")) +
    geom_hline(data = naive_lines,
               aes(yintercept = logHR_naive_bias, color = dose_selected, linetype = linetype)) +
    scale_linetype_manual(name = "Estimator Type", values = c("Debiased Estimator" = "solid", "Naive Estimator" = "dashed")) +
    labs(color = "Dose Selected") +
    ylab("Bias")

  naive_lines_sd <- est_by_dose %>%
    distinct(dose_selected, sd_logHR_naive) %>%
    mutate(linetype = "Naive Estimator")


  # Find the minimum MSE point for each dose group
  min_mse_points <- est_by_dose %>%
    group_by(dose_selected) %>%
    slice_min(mse, with_ties = FALSE)

  naive_lines_mse <- est_by_dose %>%
    distinct(dose_selected, mse_naive) %>%
    mutate(linetype = "Naive Estimator")

  # Create the plot
  mse_plot <- ggplot(est_by_dose, aes(x = w, y = sqrt(mse), color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator")) +
    geom_hline(data = naive_lines_mse,
               aes(yintercept = sqrt(mse_naive), color = dose_selected, linetype = linetype)) +
    geom_point(data = min_mse_points, aes(x = w, y = sqrt(mse)), shape = 21, fill = "white", size = 3, stroke = 1) +
    geom_text_repel(data = min_mse_points,
                    aes(x = w, y = sqrt(mse), label = paste0("min at w=", round(w, 2))),
                    size = 3,
                    nudge_y = 0.02,   # optional small upward nudge
                    direction = "y") +
    labs(color = "Dose Selected", linetype = "") +
    ylab("RMSE") +
    theme_minimal()



  # SD plot
  # sd_plot <- ggplot(est_by_dose, aes(x = w, y = sqrt(var_mu_hat_d_center), color = dose_selected, group = dose_selected)) +
  #   geom_line(aes(linetype = "Debiased Estimator")) +
  #   geom_hline(data = naive_lines_sd,
  #              aes(yintercept = sqrt(var_mu_hat_d_naive), color = dose_selected, linetype = linetype)) +
  #   ylab("SD") +
  #   labs(color = "Dose Selected")
  sd_plot <- ggplot(est_by_dose, aes(x = w, y = sd_logHR_tilde, color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator")) +
    geom_hline(data = naive_lines_sd,
               aes(yintercept = sd_logHR_naive, color = dose_selected, linetype = linetype)) +
    ylab("SD") +
    labs(color = "Dose Selected")

  naive_lines_emp_cov <- est_by_dose %>%
    distinct(dose_selected, emp_cov_naive) %>%
    mutate(linetype = "Naive Estimator")

  emp_cov_plot <- ggplot(est_by_dose, aes(x = w, y = emp_cov, color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator")) +
    geom_hline(data = naive_lines_emp_cov,
               aes(yintercept = emp_cov_naive, color = dose_selected, linetype = linetype)) +
    geom_hline(yintercept = 1 - alpha_LCL, linetype = "twodash", color = "blue") +
    labs(color = "Dose Selected") +
    ylab("Empirical Coverage")

  # Dose selection probability plot
  # dose_prob_plot <- ggplot(dose_prob_df, aes(x = Dose, y = Probability)) +
  #   geom_bar(stat = "identity", fill = "steelblue") +
  #   geom_text(aes(label = scales::percent(Probability, accuracy = 0.1)),
  #             vjust = -0.5, size = 4.5) +
  #   ylim(0, max(df$Probability) + 0.05) +
  #   ylim(0, 1) +
  #   labs(x = "Dose", y = "Selection Probability") +
  #   theme_minimal(base_size = 13)
  #
  # Combine plots

  # Extract legend from bias_plot
  legend <- cowplot::get_legend(bias_plot + theme(legend.position = "right"))

  # Remove legends from all plots
  bias_plot_clean <- bias_plot + theme(legend.position = "none")
  sd_plot_clean <- sd_plot + theme(legend.position = "none")
  mse_plot_clean <- mse_plot + theme(legend.position = "none")
  emp_cov_plot_clean <- emp_cov_plot + theme(legend.position = "none")

  # Stack plots vertically
  plots_combined <- bias_plot_clean / sd_plot_clean / mse_plot_clean / emp_cov_plot_clean +
    plot_annotation(title = title_mu)

  # Combine with legend on the right
  final_plot <- plot_grid(plots_combined, legend, ncol = 2, rel_widths = c(4, 1))

  return(final_plot)
}
