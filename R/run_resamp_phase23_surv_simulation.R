#' Runs simulations for the survival/binary resampling inference for seamless phase 2/3 trials
#'
#' @param iter
#' @param n_per_group
#' @param lambda_list
#' @param rho
#' @param p_list
#' @param alternative
#' @param alpha
#' @param alpha_LCL
#' @param ws
#' @param B
#' @param seed
#'
#' @returns
#' @export
#'
#' @examples
run_resamp_phase23_survival_simulation <- function(iter = 1000,
                                                   n_per_group_s1 = 100,
                                                   n_per_group_s2 = 100,
                                                   lambda_list = list(control = 0.1, dose1 = 0.08, dose2 = 0.06),
                                                   rho = -0.9,
                                                   p_list = list(control = 0.3, dose1 = 0.35, dose2 = 0.4),
                                                   alternative = "greater",
                                                   t0 = 12,
                                                   study_end = 36,
                                                   alpha = 0.15, alpha_LCL = 0.025,
                                                   ws = seq(0, 1, length.out = 100), B = 100,
                                                   n_iter_dose1_done = FALSE,
                                                   n_iter_dose2_done = FALSE,
                                                   seed = 123) {
  set.seed(seed)

  log_HR_tilde_all_ests <- data.frame(matrix(nrow = 0, ncol = 0))
  n_iter_dose1 = 0
  n_iter_dose2 = 0
  done = FALSE
  while (!done) {
    n_iter_dose1_done <- if(!n_iter_dose1_done) n_iter_dose1 > iter else n_iter_dose1_done
    n_iter_dose2_done <- if(!n_iter_dose2_done) n_iter_dose2 > iter else n_iter_dose2_done
    done <- (n_iter_dose1_done & n_iter_dose2_done)
    samp_inference <- run_resamp_inference_phase23_surv(n_per_group_s1 = n_per_group_s1,
                                                        n_per_group_s2 = n_per_group_s2,
                                                        lambda_list = lambda_list,
                                                        rho = rho,
                                                        p_list = p_list,
                                                        alternative = alternative,
                                                        t0 = t0,
                                                        study_end = study_end,
                                                        alpha = alpha, alpha_LCL = alpha_LCL,
                                                        ws = ws, B = B, n_iter_dose1_done = n_iter_dose1_done, n_iter_dose2_done = n_iter_dose2_done
    )


    dose_selected <- samp_inference$dose

    if (dose_selected == 1 & !n_iter_dose1_done) {
      log_HR_true <- log(lambda_list$dose1 / lambda_list$control)
      n_iter_dose1 = n_iter_dose1 + 1
      if(n_iter_dose1 %% 10 == 0) print(paste0("Dose 1:", n_iter_dose1))
    } else if (dose_selected == 2 & !n_iter_dose2_done) {
      log_HR_true <- log(lambda_list$dose2 / lambda_list$control)
      n_iter_dose2 = n_iter_dose2 + 1
      if(n_iter_dose2 %% 10 == 0) print(paste0("Dose 2:", n_iter_dose2))
    } else {
      next
    }
    logHR_naive <- samp_inference$logHR_naive
    logHR_naive_bias <- logHR_naive - log_HR_true
    var_logHR_naive <- unname(samp_inference$var_logHR_naive)
    mse_naive <- logHR_naive_bias^2
    samp_inference$logHR_w_mat <- samp_inference$logHR_w_mat %>%
      dplyr::mutate(bias = logHR_tilde - log_HR_true, mse = (logHR_tilde - log_HR_true)^2)
    logHR_s2 <- samp_inference$logHR_s2
    info_fraction <- samp_inference$info_fraction


    log_HR_tilde_all_ests <- rbind(
      log_HR_tilde_all_ests,
      cbind(samp_inference$logHR_w_mat, dose_selected, logHR_naive, logHR_naive_bias, var_logHR_naive, mse_naive, logHR_s2, info_fraction)
    )

  }
  est_by_dose <- log_HR_tilde_all_ests %>%
    dplyr::group_by(w, dose_selected) %>%
    dplyr::summarise(
      logHR_tilde_center = mean(logHR_tilde),
      sd_logHR_tilde = sd(logHR_tilde),
      bias = mean(bias),
      logHR_naive = mean(logHR_naive),
      sd_logHR_naive = mean(sqrt(var_logHR_naive)),
      logHR_naive_bias = mean(logHR_naive_bias),
      mse = mean(mse),
      mse_naive = mean(mse_naive),
      emp_cov = mean(LCL_logHR_tilde <= log_HR_true),
      emp_cov_naive = mean(LCL_logHR_naive <= log_HR_true),
      avg_info_fraction = mean(info_fraction),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dose_selected = as.character(dose_selected))

  min_w_mse <- est_by_dose %>% group_by(dose_selected) %>%
    summarise(w = w[which.min(mse)], min.mse = min(mse))
  return(list(est_by_dose = est_by_dose, min_w_mse = min_w_mse, log_HR_tilde_all_ests = log_HR_tilde_all_ests))
}
