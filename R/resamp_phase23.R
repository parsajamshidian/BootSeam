#' Function to perform resampling inference for Phase2/3 Seamless Trials
#'
#' @param mu0
#' @param mu1
#' @param mu2
#' @param sigma
#' @param n1
#' @param n2
#' @param alternative
#' @param alpha
#' @param ws
#' @param B1
#' @param B2
#' @param alpha_LCL
#' @param htest_method
#'
#' @returns
#' @export
run_resamp_inference_phase23 <- function(mu0, mu1, mu2, sigma, n1, n2, alternative,
                                         alpha = 0.15, ws = seq(0, 1, length.out = 100), B1 = 100, B2 = 100,
                                         alpha_LCL = 0.025,
                                         n_iter_dose1_done = FALSE, n_iter_dose2_done = FALSE,
                                         htest_method = "mean", input_data_list = NULL) {

  n = n1 + n2
  # Stage 1 data
  if(is.null(input_data_list)){
    y0_s1 <- rnorm(n1, mean = mu0, sd = sigma)
    y1_s1 <- rnorm(n1, mean = mu1, sd = sigma)
    y2_s1 <- rnorm(n1, mean = mu2, sd = sigma)

    y0_s2 <- rnorm(n2, mean = mu0, sd = sigma)
    y1_s2 <- rnorm(n2, mean = mu1, sd = sigma)
    y2_s2 <- rnorm(n2, mean = mu2, sd = sigma)

    dose_selection <- dose_select(y0_s1, y1_s1, y2_s1, alternative, alpha = alpha, htest_method = htest_method)
    d <- dose_selection$d



  } else{
    y0_s1 <- input_data_list$y0_s1
    y1_s1 <- input_data_list$y1_s1
    y2_s1 <- input_data_list$y2_s1

    y0_s2 <- input_data_list$y0_s2
    y1_s2 <- input_data_list$y1_s2
    y2_s2 <- input_data_list$y2_s2

    d <- input_data_list$d


  }


  if (d == 0) {
    warning("No dose was selected. Stop for futility")
    return(list(dose = d, bias_vals = NA, mu_hat_d = NA, var_mu_hat_d = NA,
                y_bar_d_biased = NA, var_y_bar_d_biased = NA, var_y_bar_d_s2 = NA))
  } else if(d == 1 & n_iter_dose1_done){
    return(list(dose = 1))
  } else if(d ==2 & n_iter_dose2_done){
    return(list(dose = 2))
  }

  if(d == 1){
    y_d_s1 = y1_s1
    y_d_s2 = y1_s2
    y_dns_s1 = y2_s1 # dose not selected
  } else{
    y_d_s1 = y2_s1
    y_d_s2 = y2_s2
    y_dns_s1 = y1_s1
  }




  data_list <- list(y0_s1 = y0_s1, y1_s1 = y1_s1, y2_s1 = y2_s1, y_d_s1 = y_d_s1, y_dns_s1 = y_dns_s1,
                    y0_s2 = y0_s2, y1_s2 = y1_s2, y2_s2 = y2_s2, y_d_s2 = y_d_s2,
                    d = d)



  # Estimates
  mu_hat_0 <- mean(c(y0_s1, y0_s2))
  var_mu_hat_0 <- var(c(y0_s1, y0_s2)) / n
  y_bar_d_s1 <- ifelse(d == 1, mean(y1_s1), mean(y2_s1))
  y_bar_0_s1 <- mean(y0_s1)
  y_bar_d_s2 <- mean(y_d_s2)
  y_bar_0_s2 <- mean(y0_s2)
  nd <- ifelse(d == 1, n1, n2)
  var_y_bar_d_s1 <- var(y_d_s1) / nd
  var_y_bar_0_s1 <- var(y0_s1) / n1
  var_y_bar_d_s2 <- var(y_d_s2) / nd
  var_y_bar_0_s2 <- var(y0_s2) / n2

  delta_hat_naive_s1 <- y_bar_d_s1 - y_bar_0_s1
  y_bar_d_naive <- mean(c(y_d_s1, y_d_s2))
  delta_hat_naive <- y_bar_d_naive - mu_hat_0
  var_delta_hat_naive = var(c(y_d_s1, y_d_s2)) / n + var_mu_hat_0
  LCL_delta_hat_naive = delta_hat_naive - qt(1 - alpha_LCL / 2, df = n - 2) * sqrt(var_delta_hat_naive)
  UCL_delta_hat_naive = delta_hat_naive + qt(1 - alpha_LCL / 2, df = n - 2) * sqrt(var_delta_hat_naive)
  # Bootstrap bias estimation
  # Bootstrap bias estimation
  bias_vals <- c()
  delta_hat_star_bs <- c()
  b1 = 0
  while (b1 <= B1) {
    boot_ind0 <- sample(1:n1, size = n1, replace = TRUE)
    boot_ind1 <- sample(1:n1, size = n1, replace = TRUE)
    boot_ind2 <- sample(1:n1, size = n1, replace = TRUE)

    y0_boot <- y0_s1[boot_ind0]
    y1_boot <- y1_s1[boot_ind1]
    y2_boot <- y2_s1[boot_ind2]

    d_boot_selection <- dose_select(y0_boot, y1_boot, y2_boot, alternative, alpha = alpha, htest_method = htest_method)
    d_boot <- d_boot_selection$d
    if (d_boot != d){
      next
    } else{
      delta_hat_nv_star_s1 <- ifelse(d_boot == 1, mean(y1_boot) - mean(y0_boot),
                                     mean(y2_boot) - mean(y0_boot))
      # bias_star = c()
      # b2 = 0
      # while(b2 <= B2){
      #   y0_boot_b2 <- y0_boot[sample(1:n1, size = n1, replace = TRUE)]
      #   y1_boot_b2 <- y1_boot[sample(1:n1, size = n1, replace = TRUE)]
      #   y2_boot_b2 <- y2_boot[sample(1:n1, size = n1, replace = TRUE)]
      #   d_boot_selection2 <- dose_select(y0_boot_b2, y1_boot_b2, y2_boot_b2, alternative, alpha = alpha, htest_method = htest_method)
      #   d_boot2 <- d_boot_selection2$d
      #   if(d_boot2 != d){
      #     next
      #   } else{
      #     delta_hat_nv_star_star_s1 <- ifelse(d_boot2 == 1, mean(y1_boot_b2) - mean(y0_boot_b2),
      #                                         mean(y2_boot_b2) - mean(y0_boot_b2))
      #     b2 = b2 + 1
      #   }
      #   bias_star = c(bias_star, delta_hat_nv_star_star_s1 - delta_hat_nv_star_s1)
      # }
      # bias_hat_star = mean(bias_star)
      #bias_vals <- c(bias_vals, delta_hat_nv_star_s1 - delta_hat_naive_s1)
      #delta_hat_star_bs = c(delta_hat_star_bs, delta_hat_nv_star_s1 - bias_hat_star)
      bias_vals <- c(bias_vals, delta_hat_nv_star_s1 - delta_hat_naive_s1)
      b1 = b1 + 1

    }
  }


  bias_hat <- mean(bias_vals)
  delta_hat_s1 <- delta_hat_naive_s1 - bias_hat # add the bias from the bootstrap
  delta_hat_s2 <- (y_bar_d_s2 - y_bar_0_s2)
  # Variance estimation
  #sigma2_1_hat <- var(delta_hat_star_bs) #+ var_y_bar_d_s1 + var_y_bar_0_s1
  sigma2_1_hat <- var(bias_vals)
  sigma2_2_hat <- var_y_bar_d_s2 + var_y_bar_0_s2


  # Create a df to store w vals
  delta_wmat <- data.frame(w = rep(NA, length(ws)),
                           delta_hat = rep(NA, length(ws)),
                           se_delta_hat = rep(NA, length(ws)),
                           LCL_T_eff = rep(NA, length(ws)),
                           UCL_T_eff = rep(NA, length(ws)))


  for(i in 1:length(ws)){
    w <- ws[i]
    delta_hat <- w * delta_hat_s1 + (1 - w) * delta_hat_s2
    var_delta_hat <- w^2 * sigma2_1_hat + (1 - w)^2 * sigma2_2_hat
    se_delta_hat <- sqrt(var_delta_hat)
    LCL_T_eff <- delta_hat - qnorm(1 - alpha_LCL/2) * se_delta_hat
    UCL_T_eff <- delta_hat + qnorm(1 - alpha_LCL/2) * se_delta_hat
    # LCL_T_eff <- delta_hat - e_u
    # UCL_T_eff <- delta_hat - e_l
    delta_wmat[i,] = c(w, delta_hat, se_delta_hat, LCL_T_eff, UCL_T_eff) #LCL_T_eff, LCL_test)
  }






  return(list(
    dose = d,
    bias_hat = bias_hat,
    delta_hat_naive = delta_hat_naive,
    var_delta_hat_naive = var_delta_hat_naive,
    LCL_delta_hat_naive = LCL_delta_hat_naive,
    UCL_delta_hat_naive = UCL_delta_hat_naive,
    delta_wmat = delta_wmat,
    data_list = data_list
  ))
}
