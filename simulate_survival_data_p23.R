# simulate_seamless_survival_data_p23 <- function(
#     n_per_group_s1 = 100,
#     n_per_group_s2 = 100,
#     lambda_list = list(control = 0.1, dose1 = 0.08, dose2 = 0.06),
#     p_list = list(p0 = 0.3, p1 = 0.4, p2 = 0.5),
#     rho = -0.9,
#     t0 = 12,
#     study_end = 36,
#     alpha = 0.15,
#     alternative = "greater"
# ) {
#   library(MASS)
#   library(dplyr)
#   # --- Stage 1 ---
#   n_stage1 <- n_per_group_s1 * 3
#   treatment_stage1 <- rep(c("control", "dose1", "dose2"), each = n_per_group_s1)
#
#   progression_probs <- ifelse(treatment_stage1 == "control", p_list$p0,
#                               ifelse(treatment_stage1 == "dose1", p_list$p1, p_list$p2))
#   mu <- c(0, 0)
#   Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
#   X <- mvrnorm(n_stage1, mu = mu, Sigma = Sigma)
#   X1 <- X[, 1]  # for progression
#   X2 <- X[, 2]  # for survival
#
#   progression <- as.integer(pnorm(X1) <= progression_probs)
#
#   lambda_stage1 <- ifelse(treatment_stage1 == "control", lambda_list$control,
#                           ifelse(treatment_stage1 == "dose1", lambda_list$dose1, lambda_list$dose2))
#   survival_time_stage1 <- qexp(pnorm(X2), rate = lambda_stage1)
#
#   entry_stage1 <- runif(n_stage1, min = 0, max = t0)
#   exit_stage1 <- entry_stage1 + survival_time_stage1
#   censoring_stage1 <- entry_stage1 + runif(n_stage1, 0, study_end - entry_stage1)
#   observed_stage1 <- pmin(exit_stage1, censoring_stage1)
#   event_stage1 <- as.integer(exit_stage1 <= censoring_stage1)
#
#   stage1_data <- data.frame(
#     treatment = treatment_stage1,
#     stage = "Stage1",
#     entry_time = entry_stage1,
#     exit_time = observed_stage1,
#     event = event_stage1,
#     progression = progression
#   )
#
#   # --- Dose Selection ---
#   y0_s1 <- stage1_data[stage1_data$treatment == "control", "progression"]
#   y1_s1 <- stage1_data[stage1_data$treatment == "dose1", "progression"]
#   y2_s1 <- stage1_data[stage1_data$treatment == "dose2", "progression"]
#
#   dose_selected <- dose_select(y0_s1, y1_s1, y2_s1,
#                                alternative = alternative,
#                                alpha = alpha,
#                                htest_method = "proportion")
#   selected_dose <- dose_selected$d
#   selected_dose_name <- paste0("dose", selected_dose)
#
#   if (selected_dose == 0) {
#     return(list(selected_dose = 0))
#   }
#
#   # --- Stage 2 (fixed hazard assignment) ---
#   n_stage2 <- n_per_group_s2 * 2
#   treatment_stage2 <- rep(c("control", selected_dose_name), each = n_per_group_s2)
#
#   lambda_stage2 <- ifelse(treatment_stage2 == "control",
#                           lambda_list$control,
#                           lambda_list[[selected_dose_name]])
#
#   entry_stage2 <- runif(n_stage2, min = t0, max = study_end)
#   survival_time_stage2 <- rexp(n_stage2, rate = lambda_stage2)
#   censoring_stage2 <- entry_stage2 + runif(n_stage2, 0, study_end - entry_stage2)
#
#   exit_stage2 <- entry_stage2 + survival_time_stage2
#   observed_stage2 <- pmin(exit_stage2, censoring_stage2)
#   event_stage2 <- as.integer(exit_stage2 <= censoring_stage2)
#
#   stage2_data <- data.frame(
#     treatment = treatment_stage2,
#     stage = "Stage2",
#     entry_time = entry_stage2,
#     exit_time = observed_stage2,
#     event = event_stage2,
#     progression = NA_integer_
#   )
#
#   # --- Combine ---
#   full_data <- rbind(stage1_data, stage2_data)
#
#   return(list(
#     full_data = full_data,
#     stage1_data = stage1_data,
#     stage2_data = stage2_data,
#     selected_dose = selected_dose
#   ))
# }

simulate_seamless_survival_data_p23 <- function(
    n_per_group_s1 = 100,
    n_per_group_s2 = 100,
    lambda_list = list(control = 0.1, dose1 = 0.08, dose2 = 0.06),
    p_list = list(p0 = 0.3, p1 = 0.4, p2 = 0.5),
    rho = -0.9,
    t0 = 12,
    study_end = 36,
    followup_stage1 = NULL,
    alpha = 0.15,
    alternative = "greater"
) {
  library(MASS)
  library(dplyr)

  if (is.null(followup_stage1)) {
    followup_stage1 <- t0
  }

  # --- Stage 1 ---
  n_stage1 <- n_per_group_s1 * 3
  treatment_stage1 <- rep(c("control", "dose1", "dose2"), each = n_per_group_s1)

  progression_probs <- ifelse(treatment_stage1 == "control", p_list$p0,
                              ifelse(treatment_stage1 == "dose1", p_list$p1, p_list$p2))
  mu <- c(0, 0)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  X <- mvrnorm(n_stage1, mu = mu, Sigma = Sigma)
  X1 <- X[, 1]  # for progression
  X2 <- X[, 2]  # for survival

  progression <- as.integer(pnorm(X1) <= progression_probs)

  lambda_stage1 <- ifelse(treatment_stage1 == "control", lambda_list$control,
                          ifelse(treatment_stage1 == "dose1", lambda_list$dose1, lambda_list$dose2))
  survival_time_stage1 <- qexp(pnorm(X2), rate = lambda_stage1)

  entry_stage1 <- runif(n_stage1, min = 0, max = t0)
  exit_stage1 <- entry_stage1 + survival_time_stage1
  censoring_stage1 <- pmin(entry_stage1 + runif(n_stage1, 0, followup_stage1), t0)
  observed_stage1 <- pmin(exit_stage1, censoring_stage1)
  event_stage1 <- as.integer(exit_stage1 <= censoring_stage1)

  stage1_data <- data.frame(
    treatment = treatment_stage1,
    stage = "Stage1",
    entry_time = entry_stage1,
    exit_time = observed_stage1,
    event = event_stage1,
    progression = progression
  )

  # --- Dose Selection ---
  y0_s1 <- stage1_data[stage1_data$treatment == "control", "progression"]
  y1_s1 <- stage1_data[stage1_data$treatment == "dose1", "progression"]
  y2_s1 <- stage1_data[stage1_data$treatment == "dose2", "progression"]

  dose_selected <- dose_select(y0_s1, y1_s1, y2_s1,
                               alternative = alternative,
                               alpha = alpha,
                               htest_method = "proportion")
  selected_dose <- dose_selected$d
  selected_dose_name <- paste0("dose", selected_dose)

  if (selected_dose == 0) {
    return(list(selected_dose = 0))
  }

  # --- Stage 2 ---
  n_stage2 <- n_per_group_s2 * 2
  treatment_stage2 <- rep(c("control", selected_dose_name), each = n_per_group_s2)

  lambda_stage2 <- ifelse(treatment_stage2 == "control",
                          lambda_list$control,
                          lambda_list[[selected_dose_name]])

  entry_stage2 <- runif(n_stage2, min = t0, max = study_end)
  survival_time_stage2 <- rexp(n_stage2, rate = lambda_stage2)
  censoring_stage2 <- entry_stage2 + runif(n_stage2, 0, study_end - entry_stage2)

  exit_stage2 <- entry_stage2 + survival_time_stage2
  observed_stage2 <- pmin(exit_stage2, censoring_stage2)
  event_stage2 <- as.integer(exit_stage2 <= censoring_stage2)

  stage2_data <- data.frame(
    treatment = treatment_stage2,
    stage = "Stage2",
    entry_time = entry_stage2,
    exit_time = observed_stage2,
    event = event_stage2,
    progression = NA_integer_
  )

  # --- Combine ---
  full_data <- rbind(stage1_data, stage2_data)

  return(list(
    full_data = full_data,
    stage1_data = stage1_data,
    stage2_data = stage2_data,
    selected_dose = selected_dose
  ))
}






set.seed(15)
sim1 <- simulate_seamless_survival_data_p23(
  n_per_group_s1 = 100,
  n_per_group_s2 = 100,
  lambda_list = list(control = 0.1, dose1 = 0.08, dose2 = 0.06),
  p_list = list(p0 = 0.3, p1 = 0.5, p2 = 0.8),
  rho = 0.9,
  t0 = 12,
  study_end = 36,
  alpha = 0.15,
  alternative = "greater"
)




sim1$stage1_data %>% group_by(treatment) %>% summarise(prop_event = mean(event), mean_entry_time = mean(entry_time),
                                                    mean_exit_time = mean(exit_time),
                                                    mean_length = mean_exit_time - mean_entry_time,
                                                    prop_prog = mean(progression))




# simulate_seamless_survival_data_independent <- function(
    #     n_per_group_s1 = 100,
#     n_per_group_s2 = 100,
#     lambda_list = list(control = 0.1, dose1 = 0.08, dose2 = 0.06),
#     p_list = list(p0 = 0.3, p1 = 0.4, p2 = 0.5),
#     t0 = 12,
#     study_end = 36,
#     alpha = 0.15,
#     alternative = "greater"
# ) {
#   library(dplyr)
#
#   # --- Stage 1 ---
#   n_stage1 <- n_per_group_s1 * 3
#   treatment_stage1 <- rep(c("control", "dose1", "dose2"), each = n_per_group_s1)
#
#   progression_probs <- ifelse(treatment_stage1 == "control", p_list$p0,
#                               ifelse(treatment_stage1 == "dose1", p_list$p1, p_list$p2))
#   progression <- rbinom(n_stage1, size = 1, prob = progression_probs)
#
#   lambda_stage1 <- ifelse(treatment_stage1 == "control", lambda_list$control,
#                           ifelse(treatment_stage1 == "dose1", lambda_list$dose1, lambda_list$dose2))
#   survival_time_stage1 <- rexp(n_stage1, rate = lambda_stage1)
#
#   entry_stage1 <- runif(n_stage1, min = 0, max = t0)
#   exit_stage1 <- entry_stage1 + survival_time_stage1
#   censoring_stage1 <- entry_stage1 + runif(n_stage1, 0, study_end - entry_stage1)
#   observed_stage1 <- pmin(exit_stage1, censoring_stage1)
#   event_stage1 <- as.integer(exit_stage1 <= censoring_stage1)
#
#   stage1_data <- data.frame(
#     treatment = treatment_stage1,
#     stage = "Stage1",
#     entry_time = entry_stage1,
#     exit_time = observed_stage1,
#     event = event_stage1,
#     progression = progression
#   )
#
#   # --- Dose Selection ---
#   y0_s1 <- stage1_data[stage1_data$treatment == "control", "progression"]
#   y1_s1 <- stage1_data[stage1_data$treatment == "dose1", "progression"]
#   y2_s1 <- stage1_data[stage1_data$treatment == "dose2", "progression"]
#
#   dose_selected <- dose_select(y0_s1, y1_s1, y2_s1,
#                                alternative = alternative,
#                                alpha = alpha,
#                                htest_method = "proportion")
#   selected_dose <- dose_selected$d
#   selected_dose_name <- paste0("dose", selected_dose)
#
#   if (selected_dose == 0) {
#     return(list(selected_dose = 0))
#   }
#
#   # --- Stage 2 ---
#   n_stage2 <- n_per_group_s2 * 2
#   treatment_stage2 <- rep(c("control", selected_dose_name), each = n_per_group_s2)
#
#   lambda_stage2 <- ifelse(treatment_stage2 == "control",
#                           lambda_list$control,
#                           lambda_list[[selected_dose_name]])
#
#   entry_stage2 <- runif(n_stage2, min = t0, max = study_end)
#   survival_time_stage2 <- rexp(n_stage2, rate = lambda_stage2)
#   censoring_stage2 <- entry_stage2 + runif(n_stage2, 0, study_end - entry_stage2)
#
#   exit_stage2 <- entry_stage2 + survival_time_stage2
#   observed_stage2 <- pmin(exit_stage2, censoring_stage2)
#   event_stage2 <- as.integer(exit_stage2 <= censoring_stage2)
#
#   stage2_data <- data.frame(
#     treatment = treatment_stage2,
#     stage = "Stage2",
#     entry_time = entry_stage2,
#     exit_time = observed_stage2,
#     event = event_stage2,
#     progression = NA_integer_
#   )
#
#   # --- Combine ---
#   full_data <- rbind(stage1_data, stage2_data)
#
#   return(list(
#     full_data = full_data,
#     stage1_data = stage1_data,
#     stage2_data = stage2_data,
#     selected_dose = selected_dose
#   ))
# }

# sim2 <- simulate_seamless_survival_data_independent(
#   n_per_group_s1 = 100,
#   n_per_group_s2 = 100,
#   lambda_list = list(control = 0.1, dose1 = 0.08, dose2 = 0.06),
#   p_list = list(p0 = 0.8, p1 = 0.6, p2 = 0.3),
#   t0 = 12,
#   study_end = 36,
#   alpha = 0.15,
#   alternative = "less"
# )


# sim2$stage1_data %>% group_by(treatment) %>% summarise(prop_event = mean(event), mean_entry_time = mean(entry_time),
#                                                        mean_exit_time = mean(exit_time),
#                                                        mean_length = mean_exit_time - mean_entry_time,
#                                                        prop_prog = mean(progression))


