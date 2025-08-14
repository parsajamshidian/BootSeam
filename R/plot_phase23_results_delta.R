#' Plots the results from the delta version of resampling phase 2/3 inference
#'
#' @param sim_result
#' @param title_mu
#' @param alpha_LCL
#'
#' @returns
#' @export
plot_phase23_results_delta <- function(sim_result, title_mu = expression(mu[0] == 0 ~ "," ~ mu[1] == 0 ~ ", "~ mu[2] == 0), alpha_LCL = 0.025) {


  large_text_theme <- theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.2, "cm"),  # wider legend lines
      legend.key.height = unit(0.6, "cm"),  # taller legend keys
      plot.title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 16)
    )




  est_by_dose <- sim_result$est_by_dose
  dose_prob_df <- sim_result$dose_prob_df

  # Prepare the data for the horizontal lines
  naive_lines <- est_by_dose %>%
    distinct(dose_selected, delta_hat_naive_bias) %>%
    mutate(linetype = "Naive Estimator")

  # Bias plot
  bias_plot <- ggplot(est_by_dose, aes(x = w, y = (bias), color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator"), lwd = 1.2) +
    geom_hline(data = naive_lines, lwd = 1.1,
               aes(yintercept = (delta_hat_naive_bias), color = dose_selected, linetype = linetype)) +
    scale_linetype_manual(name = "Estimator Type", values = c("Debiased Estimator" = "solid", "Naive Estimator" = "dashed")) +
    labs(color = "Dose Selected") +
    ylab("Bias")

  naive_lines_sd <- est_by_dose %>%
    distinct(dose_selected, sd_delta_hat_naive) %>%
    mutate(linetype = "Naive Estimator")


  # Find the minimum MSE point for each dose group
  min_mse_points <- est_by_dose %>%
    group_by(dose_selected) %>%
    slice_min(mse, with_ties = FALSE)

  naive_lines_mse <- est_by_dose %>%
    distinct(dose_selected, mse_naive) %>%
    mutate(linetype = "Naive Estimator")

  naive_lines_emp_cov <- est_by_dose %>%
    distinct(dose_selected, emp_cov_naive) %>%
    mutate(linetype = "Naive Estimator")

  # Create the plot
  mse_plot <- ggplot(est_by_dose, aes(x = w, y = sqrt(mse), color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator"), lwd = 1.2) +
    geom_hline(data = naive_lines_mse, lwd = 1.1,
               aes(yintercept = sqrt(mse_naive), color = dose_selected, linetype = linetype)) +
    geom_point(data = min_mse_points, aes(x = w, y = sqrt(mse)), shape = 21, fill = "white", size = 3, stroke = 1) +
    geom_text_repel(data = min_mse_points,
                    aes(x = w, y = sqrt(mse), label = paste0("min at w=", round(w, 2))),
                    size = 5,
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
  sd_plot <- ggplot(est_by_dose, aes(x = w, y = delta_hat_sd, color = dose_selected, group = dose_selected)) +
    geom_line(aes(linetype = "Debiased Estimator"), lwd = 1.2) +
    geom_hline(data = naive_lines_sd, lwd = 1.1,
               aes(yintercept = sd_delta_hat_naive, color = dose_selected, linetype = linetype)) +
    ylab("SD") +
    labs(color = "Dose Selected")

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
  legend <- cowplot::get_legend(bias_plot + large_text_theme + theme(legend.position = "right"))


  # Remove legends from all plots
  bias_plot_clean <- bias_plot + large_text_theme + theme(legend.position = "none")
  sd_plot_clean <- sd_plot + large_text_theme + theme(legend.position = "none")
  mse_plot_clean <- mse_plot + large_text_theme + theme(legend.position = "none")
  emp_cov_plot_clean <- emp_cov_plot + large_text_theme + theme(legend.position = "none")


  # Stack plots vertically
  plots_combined <- bias_plot_clean / sd_plot_clean / mse_plot_clean +#/ emp_cov_plot_clean +
    plot_annotation(title = title_mu, theme = theme(plot.title = element_text(size = 20, face = "bold")))

  # Combine with legend on the right
  final_plot <- plot_grid(plots_combined, legend, ncol = 2, rel_widths = c(4, 1))

  return(final_plot)
}
