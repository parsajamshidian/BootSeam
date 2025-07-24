library(survival)

#' Resampling inference for Survival/Binary Data for seamless Phase 2/3 Trials
#'
#' @param n_per_group
#' @param lambda_list
#' @param rho
#' @param p_list
#' @param alternative
#' @param alpha
#' @param alpha_LCL
#' @param ws
#' @param B
#'
#' @returns
#' @export
#'
#' @examples
run_resamp_inference_phase23_surv <- function(n_per_group = c(control = 100, dose1 = 100, dose2 = 100),
                                              lambda_list = c(control = 2, dose1 = 1, dose2 = 0.5),
                                              rho = 0.9,
                                              p_list = c(control = 0.6, dose1 = 0.3, dose2 = 0.15),
                                              alternative = "less",
                                              alpha = 0.15, alpha_LCL = 0.025,
                                              ws = seq(0, 1, length.out = 100), B = 100, n_iter_dose1_done, n_iter_dose2_done
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
  #browser()
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
  bias_vals <- c()
  b = 0
  while (b <= B) {
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

