run_resamp_phase23_survival_simulation_1 <- function(iter = 1000,
                                                     n_per_group_s1 = c(control = 100, dose1 = 100, dose2 = 100),
                                                     n_per_group_s2 = c(control = 100, dose1 = 100, dose2 = 100),
                                                     lambda = 0.1,
                                                     beta_list = c(control = 0, dose1 = log(0.8), dose2 = log(0.3)),
                                                     alternative = "less",
                                                     alpha = 0.15, alpha_LCL = 0.025,
                                                     ws = seq(0, 1, length.out = 100), B = 100,
                                                     n_iter_dose1_done = FALSE,
                                                     n_iter_dose2_done = FALSE,
                                                     seed = 123) {
  set.seed(seed)


  log_HR_tilde_all_ests <- data.frame(matrix(nrow = 0, ncol = 11))
  n_iter_dose1 = 0
  n_iter_dose2 = 0

  done = FALSE
  while (!done) {
    n_iter_dose1_done <- if(!n_iter_dose1_done) n_iter_dose1 > iter else n_iter_dose1_done
    n_iter_dose2_done <- if(!n_iter_dose2_done) n_iter_dose2 > iter else n_iter_dose2_done
    done <- (n_iter_dose1_done & n_iter_dose2_done)
    samp_inference <- run_resamp_inference_phase23_surv_1(n_per_group_s1 = n_per_group_s1,
                                                          n_per_group_s2 = n_per_group_s2,
                                                          lambda = lambda,
                                                          beta_list = beta_list,
                                                          alternative = alternative,
                                                          alpha = alpha, alpha_LCL = alpha_LCL,
                                                          ws = ws, B = B, n_iter_dose1_done = n_iter_dose1_done, n_iter_dose2_done = n_iter_dose2_done
    )


    dose_selected <- samp_inference$dose

    if (dose_selected == 1 & !n_iter_dose1_done) {
      log_HR_true <- beta_list[2]
      n_iter_dose1 = n_iter_dose1 + 1
      if(n_iter_dose1 %% 10 == 0) print(paste0("Dose 1:", n_iter_dose1))
    } else if (dose_selected == 2 & !n_iter_dose2_done) {
      log_HR_true <- beta_list[3]
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

    log_HR_tilde_all_ests <- rbind(
      log_HR_tilde_all_ests,
      cbind(samp_inference$logHR_w_mat, dose_selected, logHR_naive, logHR_naive_bias, var_logHR_naive, mse_naive)
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
      .groups = "drop"
    ) %>%
    dplyr::mutate(dose_selected = as.character(dose_selected))

  min_w_mse <- est_by_dose %>% group_by(dose_selected) %>%
    summarise(w = w[which.min(mse)], min.mse = min(mse))
  return(list(est_by_dose = est_by_dose, min_w_mse = min_w_mse, log_HR_tilde_all_ests = log_HR_tilde_all_ests))
}


run_resamp_inference_phase23_surv_1 <- function(n_per_group_s1 = c(control = 100, dose1 = 100, dose2 = 100),
                                                n_per_group_s2 = c(control = 100, dose1 = 100, dose2 = 100),
                                                lambda = 0.1,
                                                beta_list = c(control = 0, dose1 = log(0.8), dose2 = log(0.3)),
                                                alternative = "less",
                                                alpha = 0.15, alpha_LCL = 0.025,
                                                ws = seq(0, 1, length.out = 100), B = 100, n_iter_dose1_done, n_iter_dose2_done
){
  # Stage 1: simulate all groups
  stage1_data <- simulate_survival_data_general(n_per_group = n_per_group_s1,
                                                lambda = lambda,
                                                beta_list = beta_list,
                                                censor_max = 10,
                                                selected_dose = NULL)

  y0_s1 <- stage1_data[stage1_data$dose == "control", c("status")]
  y1_s1 <- stage1_data[stage1_data$dose == "dose1", c("status")]
  y2_s1 <- stage1_data[stage1_data$dose == "dose2", c("status")]

  # Estimate the log hazard ratio with a cox model
  model1 <- coxph(Surv(time, status) ~ dose, data = stage1_data)
  # model_dose1 <- coxph(Surv(time, status) ~ group, data = stage1_data %>% filter(!(group == "dose2")))
  # model_dose2 <- coxph(Surv(time, status) ~ group, data = stage1_data %>% filter(!(group == "dose1")))
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

  logHR_s1 <- ifelse(d == 1, coef(model1)["dosedose1"], coef(model1)["dosedose2"])

  # Suppose you select "doseX" for Stage 2
  selected_dose_coef_name <- paste0("dosedose", d)
  selected_dose <- paste0("dose", d)
  # Stage 2: simulate only control and selected dose
  stage2_data <- simulate_survival_data_general(n_per_group = n_per_group_s2,
                                                lambda = lambda,
                                                beta_list = beta_list,
                                                censor_max = 10,
                                                selected_dose = selected_dose)

  # run Cox model for stage 2 data
  model2 <- coxph(Surv(time, status) ~ dose, data = stage2_data)
  logHR_s2 <- coef(model2)[selected_dose_coef_name]
  # Naive estimate
  # only take selected dose from stage 1
  stage_1_data_d <- stage1_data %>% filter(dose %in% c("control", selected_dose)) %>% droplevels()
  naive_model <- coxph(Surv(time, status) ~ dose, data = rbind(stage_1_data_d, stage2_data))

  logHR_naive <- coef(naive_model)[selected_dose_coef_name]
  var_logHR_naive <- vcov(naive_model)[selected_dose_coef_name, selected_dose_coef_name]
  # Bootstrap
  bias_vals <- c()
  b = 0
  while (b <= B) {

    # Stratified bootstrap sampling
    stage1_data_boot <- stage1_data %>%
      group_by(dose) %>%
      group_modify(~ slice_sample(.x, n = n_per_group_s1[.y$dose], replace = TRUE)) %>%
      ungroup()

    # Extract status for each group
    y0_boot <- stage1_data_boot %>% filter(dose == "control") %>% pull(status)
    y1_boot <- stage1_data_boot %>% filter(dose == "dose1") %>% pull(status)
    y2_boot <- stage1_data_boot %>% filter(dose == "dose2") %>% pull(status)


    d_boot_selection <- dose_select(y0_boot, y1_boot, y2_boot, alternative, alpha = alpha, htest_method =  "proportion")
    d_boot <- d_boot_selection$d
    selected_dose_star_coef_name <- paste0("dosedose", d)
    if (d_boot != d){
      next
    } else{
      #model_boot <- coxph(Surv(time, status) ~ group, data = stage1_data_boot %>% filter(group %in% c("control", selected_dose)) %>% droplevels())
      model_boot <- coxph(Surv(time, status) ~ dose, data = stage1_data_boot)
      log_HR_star_s1 = coef(model_boot)[selected_dose_star_coef_name]
      bias_vals <- c(bias_vals, log_HR_star_s1 - logHR_s1)
      b = b + 1
    }
  }

  bias_hat <- mean(bias_vals)
  logHR_tilde_s1 <- logHR_s1 - bias_hat # add the bias from the bootstrap

  sigma2_1_hat <- var(bias_vals)
  sigma2_2_hat <- vcov(model2)[selected_dose_coef_name, selected_dose_coef_name]

  logHR_w_mat <- data.frame(w = rep(NA, length(ws)),
                            logHR_tilde = rep(NA, length(ws)),
                            var_log_HR_tilde = rep(NA, length(ws)),
                            HR_tilde = rep(NA, length(ws)),
                            LCL_logHR_tilde = rep(NA, length(ws)),
                            LCL_logHR_naive = rep(NA, length(ws)))
  LCL_logHR_naive <- logHR_naive - qnorm(1 - alpha_LCL) * sqrt(var_logHR_naive)
  for(i in 1:length(ws)){
    w <- ws[i]
    logHR_tilde <- w * logHR_tilde_s1 + (1 - w) * logHR_s2
    HR_tilde <- exp(logHR_tilde)
    var_logHR_tilde <- w^2 * sigma2_1_hat + (1 - w)^2 * sigma2_2_hat
    LCL_logHR_tilde <- logHR_tilde - qnorm(1 - alpha_LCL) * sqrt(var_logHR_tilde)
    logHR_w_mat[i,] = c(w, logHR_tilde, var_logHR_tilde, HR_tilde, LCL_logHR_tilde, LCL_logHR_naive)
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


simulate_survival_data_general <- function(n_per_group,
                                           lambda = 0.1,
                                           beta_list = NULL,
                                           censor_max = 10,
                                           selected_dose = NULL) {
  # n_per_group: named vector of sample sizes per group (e.g., c(control = 100, dose1 = 100, dose2 = 100))
  # beta_list: named vector of log hazard ratios vs control (e.g., c(dose1 = log(0.8), dose2 = log(0.6)))
  # censor_max: maximum censoring time
  # selected_dose: if specified, only simulate/control for control and selected dose (Stage 2)

  # Validate input
  if (!"control" %in% names(n_per_group)) {
    stop("You must include a 'control' group in n_per_group.")
  }

  # If selected_dose is specified, filter n_per_group and beta_list
  if (!is.null(selected_dose)) {
    if (!(selected_dose %in% names(n_per_group))) {
      stop("Selected dose must be one of the names in n_per_group.")
    }
    n_per_group <- n_per_group[c("control", selected_dose)]
    if (!is.null(beta_list)) {
      beta_list <- beta_list[selected_dose]
    }
  }

  # Create dose group vector
  dose <- unlist(mapply(rep, names(n_per_group), n_per_group))

  # Assign beta values (0 for control, others from beta_list)
  beta <- sapply(dose, function(d) {
    if (d == "control") 0 else beta_list[[d]]
  })

  # Simulate survival times
  U <- runif(length(dose))
  T <- -log(U) / (lambda * exp(beta))

  # Simulate censoring
  C <- runif(length(dose), 0, censor_max)
  time <- pmin(T, C)
  status <- as.numeric(T <= C)

  # Return data frame
  data.frame(time = time, status = status, dose = factor(dose, levels = names(n_per_group)))
}



surv_1_results <- run_resamp_phase23_survival_simulation_1(iter = 1000,
                                                           n_per_group_s1 = c(control = 100, dose1 = 100, dose2 = 100),
                                                           n_per_group_s2 = c(control = 100, dose1 = 100, dose2 = 100),
                                                           lambda = 0.1,
                                                           beta_list = c(control = 0, dose1 = log(0.9), dose2 = log(0.7)),
                                                           alternative = "less",
                                                           alpha = 0.15, alpha_LCL = 0.025,
                                                           ws = seq(0, 1, length.out = 100), B = 100,
                                                           n_iter_dose1_done = FALSE,
                                                           n_iter_dose2_done = FALSE,
                                                           seed = 100)


plot_phase23_results_surv(surv_1_results, title_mu = expression(lambda[1] == 0.09 ~ ", "~ lambda[2] == 0.07))



