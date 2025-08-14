run_resamp_inference_phase23_surv_1 <- function(n_per_group = c(control = 100, dose1 = 100, dose2 = 100),
                                              lambda_list = c(control = 2, dose1 = 1, dose2 = 0.5),
                                              rho = 0.9,
                                              p_list = c(control = 0.6, dose1 = 0.3, dose2 = 0.15),
                                              alternative = "less",
                                              alpha = 0.15, alpha_LCL = 0.025,
                                              ws = seq(0, 1, length.out = 100), B1 = 100, B2 = 100,
                                              n_iter_dose1_done = F, n_iter_dose2_done = F
){
  # Stage 1: simulate all groups
  stage1_data <- simulate_survival_data(n_per_group,
                                        lambda_list,
                                        rho,
                                        p_list
  )
  y0_s1 <- stage1_data[stage1_data$group == "control", c("status")]
  y1_s1 <- stage1_data[stage1_data$group == "dose1", c("status")]
  y2_s1 <- stage1_data[stage1_data$group == "dose2", c("status")]

  # Estimate the log hazard ratio with a cox model
  model1 <- coxph(Surv(time, status) ~ group, data = stage1_data)

  # Select dose based on proportions of survival events
  dose_selection <- dose_select(y0_s1, y1_s1, y2_s1, alternative, alpha = alpha,
                                htest_method = "proportion")
  d <- dose_selection$d

  if (d == 0) {
    warning("No dose was selected. Stop for futility")
    return(list(dose = 0))
  } else if(d == 1 & n_iter_dose1_done){
    return(list(dose = 1))
  } else if(d ==2 & n_iter_dose2_done){
    return(list(dose = 2))
  }

  logHR_s1 <- ifelse(d == 1, coef(model1)["groupdose1"], coef(model1)["groupdose2"])

  # Suppose you select "doseX" for Stage 2
  selected_dose_coef_name <- paste0("groupdose", d)
  selected_dose <- paste0("dose", d)
  # n_stage2 <- c(n2, n2)
  # names(n_stage2) <- c("control", selected_dose)
  # Stage 2: simulate only control and selected dose
  stage2_data <- simulate_survival_data(n_per_group = n_per_group[names(n_per_group) %in% c("control", selected_dose)],
                                        lambda_list,
                                        rho,
                                        p_list)

  # run Cox model for stage 2 data
  model2 <- coxph(Surv(time, status) ~ group, data = stage2_data)
  logHR_s2 <- coef(model2)[selected_dose_coef_name]
  # Naive estimate
  # only take selected dose from stage 1
  stage_1_data_d <- stage1_data %>% filter(group %in% c("control", selected_dose)) %>% droplevels()
  naive_model <- coxph(Surv(time, status) ~ group, data = rbind(stage_1_data_d, stage2_data))

  logHR_naive <- coef(naive_model)[selected_dose_coef_name]
  var_logHR_naive <- vcov(naive_model)[selected_dose_coef_name, selected_dose_coef_name]
  # Bootstrap
  bias_vals_final_est <- c()
  log_HR_hat_1_deb_bs <-  numeric(B1)
  log_HR_hat_2_bs <- numeric(B2)
  b1 = 0
  while (b1 <= B1) {
    # resample the survival data, stratified sampling
    n1 <- n_per_group[1]
    stage1_data_boot <- stage1_data %>% group_by(group) %>%
      sample_n(n1, replace = T)

    y0_boot <- stage1_data_boot[stage1_data_boot$group == "control", c("status")]
    y1_boot <- stage1_data_boot[stage1_data_boot$group == "dose1", c("status")]
    y2_boot <- stage1_data_boot[stage1_data_boot$group == "dose2", c("status")]


    d_boot_selection <- dose_select(y0_boot$status, y1_boot$status, y2_boot$status, alternative, alpha = alpha, htest_method =  "proportion")
    d_boot <- d_boot_selection$d
    selected_dose_star_coef_name <- paste0("groupdose", d)
    if (d_boot != d){
      next
    } else{
      model_boot <- coxph(Surv(time, status) ~ group, data = stage1_data_boot)
      log_HR_star_s1_b = coef(model_boot)[selected_dose_star_coef_name]
      bias_vals_final_est <- c(bias_vals_final_est, log_HR_star_s1_b - logHR_s1)
    }
    log_HR_star_star <- numeric(B2)
    b2 = 0
    while (b2 <= B2){

      stage1_data_boot2 <- stage1_data_boot %>% group_by(group) %>%
        sample_n(n1, replace = T)
      y0_boot_b2 <- stage1_data_boot2[stage1_data_boot2$group == "control", c("status")]
      y1_boot_b2 <- stage1_data_boot2[stage1_data_boot2$group == "dose1", c("status")]
      y2_boot_b2 <- stage1_data_boot2[stage1_data_boot2$group == "dose2", c("status")]
      d_boot_selection2 <- dose_select(y0_boot_b2$status, y1_boot_b2$status, y2_boot_b2$status, alternative, alpha = alpha, htest_method = "proportion")
      d_boot2 <- d_boot_selection2$d
      if(d_boot2 != d){
        next
      } else{
        model_boot2 <- coxph(Surv(time, status) ~ group, data = stage1_data_boot2)
        log_HR_star_star[b2] <- coef(model_boot2)[selected_dose_star_coef_name]
      }
      b2 = b2 + 1
    }



    bias_hat <- mean(log_HR_star_star) - log_HR_star_s1_b
    log_HR_hat_1_deb_bs[b1] <- log_HR_star_s1_b - bias_hat

    # Second stage estimate
    stage2_data_boot <- stage2_data %>% group_by(group) %>%
      sample_n(n2, replace = T)

    y0_boot_s2 <-  stage2_data_boot[stage2_data_boot$group == "control", c("status")]
    yd_boot_s2 <-  stage2_data_boot[stage2_data_boot$group == selected_dose, c("status")]

    model_boot_s2 <- coxph(Surv(time, status) ~ group, data = stage2_data_boot)
    log_HR_hat_2_bs[b1] <- coef(model_boot_s2)[selected_dose_star_coef_name]

    b1 = b1 + 1


  }

  logHR_tilde_s1 <- logHR_s1 - mean(bias_vals_final_est) # add the bias from the bootstrap

  sigma2_1_hat <- var(bias_vals_final_est)
  sigma2_2_hat <- vcov(model2)[selected_dose_coef_name, selected_dose_coef_name]

  empty_ws_vec = rep(NA, length(ws))
  logHR_w_mat <- data.frame(w = empty_ws_vec,
                            logHR_tilde_final_est = empty_ws_vec,
                            var_logHR_tilde = empty_ws_vec,
                            LCL = empty_ws_vec,
                            UCL = empty_ws_vec,
                            LCL_logHR_naive = rep(NA, length(ws)),
                            UCL_logHR_naive = rep(NA, length(ws)))
  LCL_logHR_naive <- logHR_naive - qnorm(1 - alpha_LCL) * sqrt(var_logHR_naive)
  UCL_logHR_naive <- logHR_naive + qnorm(1 - alpha_LCL) * sqrt(var_logHR_naive)
  for(i in 1:length(ws)){
    w <- ws[i]
    logHR_tilde_bs <- w * log_HR_hat_1_deb_bs + (1 - w) * log_HR_hat_2_bs
    logHR_tilde_final_est <- w * logHR_tilde_s1 + (1 - w) * logHR_s2
    var_logHR_tilde <- w^2 * sigma2_1_hat + (1 - w)^2 * sigma2_2_hat
    LCL <- quantile(logHR_tilde_bs, probs = alpha_LCL / 2)
    UCL <- quantile(logHR_tilde_bs, probs = 1 - alpha_LCL/2)

    logHR_w_mat[i,] = c(w, logHR_tilde_final_est, var_logHR_tilde, LCL, UCL, LCL_logHR_naive, UCL_logHR_naive)
  }

  return(list(
    dose = d,
    bias_hat = bias_hat,
    logHR_naive = logHR_naive,
    var_logHR_naive = var_logHR_naive,
    logHR_tilde_s1 = logHR_tilde_s1,
    logHR_s2 = logHR_s2,
    logHR_w_mat = logHR_w_mat
  ))





}

run_resamp_phase23_survival_simulation_new <- function(iter = 1000,
                                                   n_per_group = c(control = 100, dose1 = 100, dose2 = 100),
                                                   lambda_list = c(control = 2, dose1 = 1, dose2 = 0.5),
                                                   rho = 0.9,
                                                   p_list = c(control = 0.6, dose1 = 0.3, dose2 = 0.15),
                                                   alternative = "less",
                                                   alpha = 0.15, alpha_LCL = 0.025,
                                                   ws = seq(0, 1, length.out = 100), B1 = 100, B2 = 100,
                                                   seed = 123) {
  set.seed(seed)

  log_HR_tilde_all_ests <- data.frame(matrix(nrow = 0, ncol = 0))
  n_iter_dose1 = 0
  n_iter_dose2 = 0
  n_iter_dose1_done = FALSE
  n_iter_dose2_done = FALSE
  done = FALSE
  for(i in 1:iter){
  # while (!done) {
  #   n_iter_dose1_done <- n_iter_dose1 > iter
  #   n_iter_dose2_done <- n_iter_dose2 > iter
  #   done <- (n_iter_dose1_done & n_iter_dose2_done)
    samp_inference <- run_resamp_inference_phase23_surv_1(n_per_group = n_per_group,
                                                        lambda_list = lambda_list,
                                                        rho = rho,
                                                        p_list = p_list,
                                                        alternative = alternative,
                                                        alpha = alpha, alpha_LCL = alpha_LCL,
                                                        ws = ws, B1 = B1, B2 = B2,
                                                        n_iter_dose1_done = n_iter_dose1_done, n_iter_dose2_done = n_iter_dose2_done
    )


    dose_selected <- samp_inference$dose

    if (dose_selected == 1 & !n_iter_dose1_done) {
      log_HR_true <- log(lambda_list["dose1"] / lambda_list[1])
      n_iter_dose1 = n_iter_dose1 + 1
      if(n_iter_dose1 %% 100 == 0) print(paste0("Dose 1:", n_iter_dose1))
    } else if (dose_selected == 2 & !n_iter_dose2_done) {
      log_HR_true <- log(lambda_list["dose2"] / lambda_list[1])
      n_iter_dose2 = n_iter_dose2 + 1
      if(n_iter_dose2 %% 100 == 0) print(paste0("Dose 2:", n_iter_dose2))
    } else {
      next
    }
    logHR_naive <- samp_inference$logHR_naive
    logHR_naive_bias <- logHR_naive - log_HR_true
    var_logHR_naive <- unname(samp_inference$var_logHR_naive)
    mse_naive <- logHR_naive_bias^2
    samp_inference$logHR_w_mat <- samp_inference$logHR_w_mat %>%
      dplyr::mutate(bias = logHR_tilde_final_est - log_HR_true, mse = (logHR_tilde_final_est - log_HR_true)^2)

    log_HR_tilde_all_ests <- rbind(
      log_HR_tilde_all_ests,
      cbind(samp_inference$logHR_w_mat, dose_selected, logHR_naive, logHR_naive_bias, var_logHR_naive, mse_naive)
    )
  message("Completed iteration: ", i)
  }
  est_by_dose <- log_HR_tilde_all_ests %>%
    dplyr::group_by(w, dose_selected) %>%
    dplyr::summarise(
      logHR_tilde_center = mean(logHR_tilde_final_est),
      sd_logHR_tilde = sd(logHR_tilde_final_est),
      bias = mean(bias),
      logHR_naive = mean(logHR_naive),
      sd_logHR_naive = mean(sqrt(var_logHR_naive)),
      logHR_naive_bias = mean(logHR_naive_bias),
      mse = mean(mse),
      mse_naive = mean(mse_naive),
      emp_cov = mean(LCL <= log_HR_true & UCL >= log_HR_true),
      emp_cov_naive = mean(LCL_logHR_naive <= log_HR_true & UCL_logHR_naive >= log_HR_true),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dose_selected = as.character(dose_selected))

  min_w_mse <- est_by_dose %>% group_by(dose_selected) %>%
    summarise(w = w[which.min(mse)], min.mse = min(mse))
  return(list(est_by_dose = est_by_dose, min_w_mse = min_w_mse, log_HR_tilde_all_ests = log_HR_tilde_all_ests))
}




simulate_survival_data <- function(n_per_group,
                                   lambda_list = NULL,
                                   rho = NULL,
                                   p_list = NULL) {
  # n_per_group: named vector of sample sizes per group (e.g., c(control = 100, dose1 = 100, dose2 = 100))
  # selected_dose: if specified, only simulate for control and selected dose
  df_surv <- data.frame()
  for(i in 1:(length(n_per_group))){
    n_g <- n_per_group[i]
    group_name <- names(n_per_group)[i]
    Z1 = rnorm(n_g)
    Z2 = rnorm(n_g)
    # Construct correlated variables
    X1 <- Z1
    X2 <- (rho * Z1 + sqrt(1 - rho^2) * Z2)

    Y1 <- as.numeric(pnorm(X1) <= p_list[i])
    Y2 <- qexp(pnorm(X2), rate = lambda_list[i])

    df_surv <- rbind(df_surv,
                     data.frame(
                       group = group_name,
                       status = Y1,
                       time = Y2
                     ))

  }
  return(df_surv)
}



surv_sim_results1 <- run_resamp_phase23_survival_simulation_new(iter = 10,
                                                           n_per_group = c(control = 100, dose1 = 100, dose2 = 100),
                                                           lambda_list = c(control = 0.1, dose1 = 0.08, dose2 = 0.05),
                                                           rho = 0.9,
                                                           p_list = c(control = 0.5, dose1 = 0.4, dose2 = 0.3),
                                                           alternative = "less",
                                                           alpha = 0.15, alpha_LCL = 0.05,
                                                           ws = seq(0, 1, length.out = 100), B1 = 100, B2 = 50,
                                                           seed = 10)


