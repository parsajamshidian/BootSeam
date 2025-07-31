run_resamp_inference_phase23_1 <- function(mu0, mu1, mu2, sigma, n1, n2, alternative,
                                           alpha = 0.15, ws = seq(0, 1, length.out = 100), B1 = 100, B2 = 100,
                                           alpha_LCL = 0.025,
                                           n_iter_dose1_done = FALSE, n_iter_dose2_done = FALSE,
                                           htest_method = "mean") {


  n = n1 + n2
  # Stage 1 data
  y0_s1 <- rnorm(n1, mean = mu0, sd = sigma)
  y1_s1 <- rnorm(n1, mean = mu1, sd = sigma)
  y2_s1 <- rnorm(n1, mean = mu2, sd = sigma)

  dose_selection <- dose_select(y0_s1, y1_s1, y2_s1, alternative, alpha = alpha, htest_method = htest_method)
  d <- dose_selection$d

  if (d == 0) {
    warning("No dose was selected. Stop for futility")
    return(list(dose = d, bias_vals = NA, mu_hat_d = NA, var_mu_hat_d = NA,
                y_bar_d_biased = NA, var_y_bar_d_biased = NA, var_y_bar_d_s2 = NA))
  } else if(d == 1 & n_iter_dose1_done){
    return(list(dose = 1))
  } else if(d ==2 & n_iter_dose2_done){
    return(list(dose = 2))
  }

  # Stage 2 data
  y0_s2 <- rnorm(n2, mean = mu0, sd = sigma)
  if(d == 1){
    y_d_s1 = y1_s1
    y_d_s2 = rnorm(n2, mean = mu1, sd = sigma)
    y_dns_s1 = y2_s1 # dose not selected
  } else{
    y_d_s1 = y2_s1
    y_d_s2 = rnorm(n2, mean = mu2, sd = sigma)
    y_dns_s1 = y1_s1
  }

  # variance estimates
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
  delta_hat_naive <- mean(y_bar_d_naive - mu_hat_0)
  var_delta_hat_naive = var(c(y_d_s1, y_d_s2)) / n + var_mu_hat_0
  LCL_delta_hat_naive = delta_hat_naive - qt(1 - alpha_LCL / 2, df = n - 2) * sqrt(var_delta_hat_naive)
  UCL_delta_hat_naive = delta_hat_naive + qt(1 - alpha_LCL / 2, df = n - 2) * sqrt(var_delta_hat_naive)
  # Bootstrap bias estimation



  delta_hat_1_deb_bs <- numeric(B1)
  delta_hat_2_bs <- numeric(B1)

  b1 = 0
  while (b1 < B1) {

    # Resample stage 1 data
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
      delta_1_star_b <- ifelse(d_boot == 1, mean(y1_boot) - mean(y0_boot),
                               mean(y2_boot) - mean(y0_boot))
    }

    # Now do nested bootstrap for bias correction
    delta_star_star <- numeric(B2)
    b2 = 0
    while (b2 < B2) {
      y0_boot_b2 <- y0_boot[sample(1:n1, size = n1, replace = TRUE)]
      y1_boot_b2 <- y1_boot[sample(1:n1, size = n1, replace = TRUE)]
      y2_boot_b2 <- y2_boot[sample(1:n1, size = n1, replace = TRUE)]
      d_boot_selection2 <- dose_select(y0_boot_b2, y1_boot_b2, y2_boot_b2, alternative, alpha = alpha, htest_method = htest_method)
      d_boot2 <- d_boot_selection2$d
      if(d_boot2 != d){
        next
      } else{
        delta_star_star[b2] <- ifelse(d_boot2 == 1, mean(y1_boot_b2) - mean(y0_boot_b2),
                                      mean(y2_boot_b2) - mean(y0_boot_b2))
      }
      b2 = b2 + 1
    }
    browser()
    # Estimate bias within this outer resample
    bias_hat <- mean(delta_star_star) - delta_1_star_b

    # Debiased estimate
    delta_hat_1_deb_bs[b1] <- delta_1_star_b - bias_hat

    # Second stage estimate
    boot_ind0_s2 <- sample(1:n2, size = n2, replace = TRUE)
    boot_indd_s2 <- sample(1:n2, size = n2, replace = TRUE)

    y0_boot_s2 <- y0_s2[boot_ind0_s2]
    yd_boot_s2 <- y_d_s2[boot_indd_s2]

    delta_hat_2_bs[b1] <- mean(yd_boot_s2) - mean(y0_boot_s2)
    b1 = b1 + 1

  }




  empty_ws_vec = rep(NA, length(ws))
  delta_ci_mat = data.frame(w = empty_ws_vec,
                            delta_hat = empty_ws_vec,
                            LCL = empty_ws_vec,
                            UCL = empty_ws_vec)
  for (i in 1:length(ws)) {
    w <- ws[i]
    delta_hat_full_bs <- w * delta_hat_1_deb_bs + (1 - w) * delta_hat_2_bs
    LCL <- quantile(delta_hat_full_bs, probs = alpha_LCL / 2)
    UCL <- quantile(delta_hat_full_bs, probs = 1 - alpha_LCL/2)
    delta_hat <- mean(delta_hat_full_bs)
    delta_ci_mat[i, ] <- c(w, delta_hat, LCL, UCL)
    # Store CI, check coverage, etc.
  }

  return(list(
    dose = d,
    delta_ci_mat = delta_ci_mat,
    delta_hat_naive = delta_hat_naive,
    LCL_delta_hat_naive = LCL_delta_hat_naive,
    UCL_delta_hat_naive = UCL_delta_hat_naive
  ))
}
