#' Runs a simulation over n iterations for the resampling phase 2/3 inference
#'
#' @param iter
#' @param mu0
#' @param mu1
#' @param mu2
#' @param n1
#' @param n2
#' @param ws
#' @param alpha
#' @param alpha_LCL
#' @param htest_method
#' @param seed
#'
#' @returns
#' @export
#'
#' @examples
run_resamp_phase23_simulation <- function(iter = 1000,
                                          mu0 = 0,
                                          mu1 = 0,
                                          mu2 = 0,
                                          n1 = 100,
                                          n2 = 100,
                                          B1 = 100,
                                          B2 = 100,
                                          B_errs = 200,
                                          ws = seq(0, 1, length.out = 100),
                                          alpha = 0.15,
                                          alpha_LCL = 0.025,
                                          htest_method = "mean",
                                          seed = 123) {
  set.seed(seed)

  delta_all_ests <- data.frame(matrix(nrow = 0, ncol = 10))

  n_iter_dose1 = 0
  n_iter_dose2 = 0
  n_iter_dose1_done = FALSE
  n_iter_dose2_done = FALSE
  done = FALSE
  while(!done) {
    n_iter_dose1_done <- n_iter_dose1 > iter
    n_iter_dose2_done <- n_iter_dose2 > iter
    done <- (n_iter_dose1_done & n_iter_dose2_done)
    samp_inference <- run_resamp_inference_phase23(
      mu0 = mu0, mu1 = mu1, mu2 = mu2,
      sigma = 1, n1 = n1, n2 = n2,
      alternative = "greater",
      alpha = alpha, ws = ws, B1 = B1, B2 = B2,
      alpha_LCL = alpha_LCL,
      n_iter_dose1_done = n_iter_dose1_done, n_iter_dose2_done = n_iter_dose2_done,
      htest_method = htest_method
    )

    dose_selected <- samp_inference$dose

    if (dose_selected == 1 & !n_iter_dose1_done) {
      delta_true <- mu1 - mu0
      n_iter_dose1 = n_iter_dose1 + 1
      if(n_iter_dose1 %% 10 == 0) print(paste0("Dose1:",n_iter_dose1))
    } else if (dose_selected == 2 & !n_iter_dose2_done) {
      delta_true <- mu2 - mu0
      n_iter_dose2 = n_iter_dose2 + 1
      if(n_iter_dose2 %% 10 == 0) print(paste0("Dose2:",n_iter_dose2))
    } else {
      next
    }

    mu_hat_0 <- samp_inference$mu_hat_0
    delta_hat_naive <- samp_inference$delta_hat_naive
    delta_hat_naive_bias <- delta_hat_naive - delta_true
    LCL_delta_hat_naive <- samp_inference$LCL_delta_hat_naive
    UCL_delta_hat_naive <- samp_inference$UCL_delta_hat_naive
    var_delta_hat_naive <- samp_inference$var_delta_hat_naive
    mse_naive <- delta_hat_naive_bias^2
    samp_inference$delta_wmat <- samp_inference$delta_wmat %>%
      dplyr::mutate(bias = delta_hat - delta_true, mse = (delta_hat - delta_true)^2)

    delta_all_ests <- rbind(
      delta_all_ests,
      cbind(samp_inference$delta_wmat, dose_selected, delta_hat_naive, delta_hat_naive_bias, var_delta_hat_naive, mse_naive, LCL_delta_hat_naive)
    )
  }
  est_by_dose <- delta_all_ests %>%
    dplyr::group_by(w, dose_selected) %>%
    dplyr::summarise(
      delta_hat_center = mean(delta_hat),
      delta_hat_sd = sd(delta_hat),
      delta_hat_se = mean(se_delta_hat),
      bias = mean(bias),
      delta_hat_naive_center = mean(delta_hat_naive),
      delta_hat_naive_bias = mean(delta_hat_naive_bias),
      var_delta_hat_naive = mean(var_delta_hat_naive),
      sd_delta_hat_naive = sd(delta_hat_naive),
      mse = mean(mse),
      mse_naive = mean(mse_naive),
      emp_cov = mean(LCL_T_eff <= delta_true & UCL_T_eff >= delta_true),
      emp_cov_naive = mean(LCL_delta_hat_naive <= delta_true & UCL_delta_hat_naive >= delta_true),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dose_selected = as.character(dose_selected))


  min_w_mse <- est_by_dose %>% group_by(dose_selected) %>%
    summarise(w = w[which.min(mse)], min.mse = min(mse))
  return(list(est_by_dose = est_by_dose, min_w_mse = min_w_mse, delta_all_ests = delta_all_ests))
}
