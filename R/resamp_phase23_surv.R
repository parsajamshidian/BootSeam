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
run_resamp_inference_phase23_surv <- function(n_per_group_s1 = 100,
                                              n_per_group_s2 = 100,
                                              lambda_list = list(control = 0.1, dose1 = 0.08, dose2 = 0.06),
                                              rho = 0.9,
                                              p_list = list(p0 = 0.5, p1 = 0.4, p2 = 0.3),
                                              alternative = "less",
                                              t0 = 12,
                                              study_end = 36,
                                              alpha = 0.15, alpha_LCL = 0.025,
                                              ws = seq(0, 1, length.out = 100), B = 100, n_iter_dose1_done = F, n_iter_dose2_done = F
){
  # Stage 1: simulate all groups
  survival_data <- simulate_seamless_survival_data_p23(
    n_per_group_s1 = n_per_group_s1,
    n_per_group_s2 = n_per_group_s2,
    lambda_list = lambda_list,
    p_list = p_list,
    rho = rho,
    t0 = t0,
    study_end = study_end,
    alpha = alpha,
    alternative = alternative
  )
  # stage1_data <- generate_survival_data(n_per_group_s1,
  #                                       lambda_list,
  #                                       rho,
  #                                       p_list
  # )
  #
  # y0_s1 <- stage1_data[stage1_data$group == "control", c("status")]
  # y1_s1 <- stage1_data[stage1_data$group == "dose1", c("status")]
  # y2_s1 <- stage1_data[stage1_data$group == "dose2", c("status")]
  #
  # # Estimate the log hazard ratio with a cox model
  # model_dose1 <- coxph(Surv(time, status) ~ group, data = stage1_data %>% filter(!(group == "dose2")))
  # model_dose2 <- coxph(Surv(time, status) ~ group, data = stage1_data %>% filter(!(group == "dose1")))
  # Select dose based on proportions of survival events
  # dose_selection <- dose_select(y0_s1, y1_s1, y2_s1, alternative, alpha = alpha,
  #                               htest_method = "proportion")
  d <- survival_data$selected_dose
  if (d == 0) {
    warning("No dose was selected. Stop for futility")
    return(list(dose = 0))
  } else if(d == 1 & n_iter_dose1_done){
    return(list(dose = 1))
  } else if(d ==2 & n_iter_dose2_done){
    return(list(dose = 2))
  }
  stage1_data = survival_data$stage1_data

  # Suppose you select "doseX" for Stage 2
  selected_dose_coef_name <- paste0("treatmentdose", d)
  selected_dose <- paste0("dose", d)
  # Stage 2: simulate only control and selected dose
  # stage2_data <- generate_survival_data(n_per_group = n_per_group_s2[names(n_per_group_s2) %in% c("control", selected_dose)],
  #                                       lambda_list,
  #                                       rho,
  #                                       p_list)
  # Stage 1 model
  model1 <- coxph(Surv(entry_time, exit_time, event) ~ treatment, data = stage1_data)

  logHR_s1 <- ifelse(d == 1, coef(model1)["treatmentdose1"], coef(model1)["treatmentdose2"])

  stage2_data <- survival_data$stage2_data
  # run Cox model for stage 2 data
  model2 <- coxph(Surv(entry_time, exit_time, event) ~ treatment, data = stage2_data)

  # Filter Stage 1 data to selected dose and control
  stage1_subset <- stage1_data %>%
    filter(treatment %in% c("control", selected_dose))

  # Count number of events in Stage 1 for selected arms
  n_events_stage1 <- sum(stage1_subset$event)

  # Filter Stage 2 data to selected dose and control
  stage2_subset <- stage2_data %>%
    filter(treatment %in% c("control", selected_dose))

  # Count number of events in Stage 2 for selected arms
  n_events_stage2 <- sum(stage2_subset$event)

  # Compute information fraction
  info_fraction <- n_events_stage1 / (n_events_stage1 + n_events_stage2)


  logHR_s2 <- coef(model2)[selected_dose_coef_name]
  # Naive estimate
  df_all <- survival_data$full_data
  naive_model <- coxph(Surv(entry_time, exit_time, event) ~ treatment, data = df_all)
  logHR_naive <- coef(naive_model)[selected_dose_coef_name]
  var_logHR_naive <- vcov(naive_model)[selected_dose_coef_name, selected_dose_coef_name]
  # Bootstrap
  bias_vals <- c()
  b = 0
  while (b <= B) {

    # Stratified bootstrap sampling
    stage1_data_boot <- stage1_data %>%
      group_by(treatment) %>%
      group_modify(~ slice_sample(.x, n = n_per_group_s1, replace = TRUE)) %>%
      ungroup()

    # Extract status for each group
    y0_boot <- stage1_data_boot %>% filter(treatment == "control") %>% pull(progression)
    y1_boot <- stage1_data_boot %>% filter(treatment == "dose1") %>% pull(progression)
    y2_boot <- stage1_data_boot %>% filter(treatment == "dose2") %>% pull(progression)


    d_boot_selection <- dose_select(y0_boot, y1_boot, y2_boot, alternative, alpha = alpha, htest_method =  "proportion")
    d_boot <- d_boot_selection$d
    selected_dose_star_coef_name <- paste0("treatmentdose", d)
    if (d_boot != d){
      next
    } else{
      #model_boot <- coxph(Surv(time, status) ~ group, data = stage1_data_boot %>% filter(group %in% c("control", selected_dose)) %>% droplevels())
      model_boot <- coxph(Surv(entry_time, exit_time, event) ~ treatment, data = stage1_data_boot)
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
    logHR_w_mat = logHR_w_mat,
    info_fraction = info_fraction
  ))





}




