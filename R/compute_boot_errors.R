#' Compute boot errors for treatment effect estimation
#'
#' @param resamp_results
#' @param n1
#' @param n2
#' @param B_errs
#' @param alpha
#' @param htest_method
#' @param alpha_LCL
#'
#' @returns
#' @export
#'
#' @examples
compute_delta_boot_error_quantiles <- function(resamp_results, n1, n2, alternative = "greater",
                                B_errs = 100, alpha = 0.15, htest_method = "mean", alpha_LCL = 0.05) {
  # Intervals
  # Extract the data df from the original run
  data_list <- resamp_results$data_list
  d <- resamp_results$dose

  # Combine the simulated data to form the "true" distribution
  Y_0 = c(data_list$y0_s1, data_list$y0_s2)
  Y_d = c(data_list$y_d_s1, data_list$y_d_s2)
  Y_dns = data_list$y_dns_s1
  delta_star_true = mean(Y_d) - mean(Y_0)
  delta_boot_error_mat = data.frame(matrix(nrow = 0, ncol = 0))
  b1 = 0
  while(b1 < B_errs){
    Y_0_s1_boot <- Y_0[sample(1:length(Y_0), replace = TRUE)]
    Y_d_s1_boot <- Y_d[sample(1:length(Y_d), replace = TRUE)]
    Y_dns_s1_boot <- Y_dns[sample(1:length(Y_dns), replace = TRUE)]

    if(d == 1){
      dose_selection_boot <- dose_select(Y_0_s1_boot, Y_d_s1_boot, Y_dns_s1_boot, alternative = alternative, alpha = alpha, htest_method = htest_method)
      } else if (d == 2){
      dose_selection_boot <- dose_select(Y_0_s1_boot, Y_dns_s1_boot, Y_d_s1_boot, alternative = alternative, alpha = alpha, htest_method = htest_method)
    }
    d_boot <- dose_selection_boot$d
    if(d_boot != d){
      next
    } else{
      # sample stage 2 data
      Y_0_s2_boot <- Y_0[sample(1:length(Y_0), replace = TRUE)]
      Y_d_s2_boot <- Y_d[sample(1:length(Y_d), replace = TRUE)]
      if(d == 1){
        input_data_list <- list(y0_s1 = Y_0_s1_boot, y1_s1 = Y_d_s1_boot, y2_s1 = Y_dns_s1_boot,
                                y0_s2 = Y_0_s2_boot, y1_s2 = Y_d_s2_boot, y2_s2 = NULL, d = d)
      } else if(d == 2){
        input_data_list <- list(y0_s1 = Y_0_s1_boot, y1_s1 =Y_dns_s1_boot , y2_s1 = Y_d_s1_boot,
                                y0_s2 = Y_0_s2_boot, y1_s2 = NULL, y2_s2 = Y_d_s2_boot, d = d)
      }



      resamp_results_boot <- run_resamp_inference_phase23(mu0 = NULL, mu1 = NULL, mu2 = NULL, sigma = NULL, n1 = n1, n2 = n2, alternative = alternative,
                                                     alpha = alpha, ws = seq(0, 1, length.out = 100), B1 = 100, B2 = 100,
                                                     alpha_LCL = alpha_LCL,
                                                     htest_method = "mean", input_data_list = input_data_list)

      # compute the bootstrapped estimation errors

      delta_wmat_boot <- resamp_results_boot$delta_wmat

      delta_wmat_boot <- delta_wmat_boot %>% mutate(boot_errors = delta_hat - delta_star_true)

      delta_boot_error_mat = rbind(delta_boot_error_mat, delta_wmat_boot)
      b1 = b1 + 1
    }
  }
  delta_boot_qtl_mat <- delta_boot_error_mat %>% group_by(w) %>%
    summarise(e_l = quantile(boot_errors, probs = alpha_LCL / 2),
              e_u = quantile(boot_errors, probs = 1 - alpha_LCL / 2)
              )

  return(delta_boot_qtl_mat)


}

