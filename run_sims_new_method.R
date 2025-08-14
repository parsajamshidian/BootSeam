library(tidyverse)
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
  bias_vals_final_s1_est <- c()
  while (b1 <= B1) {

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
    bias_vals_final_s1_est <- c(bias_vals_final_s1_est, delta_1_star_b - delta_hat_naive_s1)

    # Now do nested bootstrap for bias correction
    delta_star_star <- numeric(B2)
    b2 = 0
    while (b2 <= B2) {
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



  delta_hat_s1 <- delta_hat_naive_s1 - mean(bias_vals_final_s1_est) # add the bias from the bootstrap
  delta_hat_s2 <- (y_bar_d_s2 - y_bar_0_s2)

  empty_ws_vec = rep(NA, length(ws))
  delta_ci_mat = data.frame(w = empty_ws_vec,
                            delta_hat = empty_ws_vec,
                            delta_hat_final_est = empty_ws_vec,
                            LCL = empty_ws_vec,
                            UCL = empty_ws_vec)
  for (i in 1:length(ws)) {
    w <- ws[i]
    delta_hat_full_bs <- w * delta_hat_1_deb_bs + (1 - w) * delta_hat_2_bs
    delta_hat_final_est <- w * delta_hat_s1 + (1 - w) * delta_hat_s2
    LCL <- quantile(delta_hat_full_bs, probs = alpha_LCL / 2)
    UCL <- quantile(delta_hat_full_bs, probs = 1 - alpha_LCL/2)
    delta_hat <- mean(delta_hat_full_bs)
    delta_ci_mat[i, ] <- c(w, delta_hat, delta_hat_final_est, LCL, UCL)
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

dose_select <- function(y0, y1, y2, alternative, alpha = 0.15, htest_method = "mean") {

  if(htest_method == "mean"){
    t1 <- t.test(y1, y0, alternative = alternative)
    t2 <- t.test(y2, y0, alternative = alternative)
  } else if (htest_method == "proportion") {
    t1 <- prop.test(c(sum(y1), sum(y0)), c(length(y1), length(y0)), alternative = alternative)
    t2 <- prop.test(c(sum(y2), sum(y0)), c(length(y2), length(y0)), alternative = alternative)
  }else if (htest_method == "Dunnett"){
    groups <- rep(c("control", "dose1", "dose2"), each = length(y0))
    Dtest <- DunnettTest(c(y0, y1, y2), groups)
    p_vals <- Dtest$control[, 4]
    t1 <- list(p.value = p_vals[1])
    t2 <- list(p.value = p_vals[2])
  } else {
    stop("Invalid htest_method specified")
  }

  fw_alpha <- alpha
  if (min(t1$p.value, t2$p.value) < fw_alpha) {
    if (t2$p.value < t1$p.value) {
      d <- 2
    } else {
      d <- 1
    }
  } else {
    d <- 0
  }
  return(list(d = d, t1_pval = t1$p.value, t2_pval = t2$p.value))
}


run_phase23_simulation_new <- function(n_sim = 100,
                                   mu0 = 0, mu1 = 0.2, mu2 = 0.3, sigma = 1,
                                   n1 = 100, n2 = 100,
                                   alpha = 0.15, alpha_LCL = 0.025,
                                   ws = seq(0, 1, length.out = 500),
                                   B1 = 100, B2 = 100,
                                   htest_method = "mean",
                                   seed = 7) {
  set.seed(seed)
  delta_all_ests <- data.frame()
  for (i in 1:n_sim) {
    results <- run_resamp_inference_phase23_1(
      mu0 = mu0, mu1 = mu1, mu2 = mu2, sigma = sigma,
      n1 = n1, n2 = n2,
      alternative = "greater",
      alpha = alpha, ws = ws,
      B1 = B1, B2 = B2, alpha_LCL = alpha_LCL,
      n_iter_dose1_done = FALSE, n_iter_dose2_done = FALSE,
      htest_method = htest_method
    )

    dose_selected <- results$dose
    if (dose_selected == 0) next
    delta_true <- ifelse(dose_selected == 1, mu1 - mu0, mu2 - mu0)
    delta_hat_naive <- results$delta_hat_naive
    LCL_naive <- results$LCL_delta_hat_naive
    UCL_naive <- results$UCL_delta_hat_naive

    covered_naive <- as.integer(LCL_naive <= delta_true & UCL_naive >= delta_true)
    naive_bias <- (delta_hat_naive - delta_true)

    delta_ci_mat <- results$delta_ci_mat %>%
      mutate(
        covered = as.integer(LCL <= delta_true & UCL >= delta_true),
        bias = (delta_hat_final_est - delta_true),
        dose_selected = dose_selected,
        delta_hat_naive = delta_hat_naive,
        naive_bias = naive_bias,
        covered_naive = covered_naive
      )

    delta_all_ests <- bind_rows(delta_all_ests, delta_ci_mat)

    message("Completed iteration: ", i)
  }

  return(delta_all_ests)
}

delta_all_ests <- run_phase23_simulation_new(n_sim = 1000, mu0 = 0, mu1 = 0, mu2 = 0, n1 = 100, n2 = 100, B1 = 500,
                                             seed = 5, alpha_LCL = 0.05, alpha = 0.15)



delta_all_ests_summarized <- delta_all_ests %>% mutate(dose_selected = as.character(dose_selected)) %>%
  group_by(w, dose_selected) %>%
  summarise(emp_cov = mean(covered),
            bias = mean(bias),
            naive_center = mean(delta_hat_naive),
            naive_bias = mean(naive_bias),
            emp_cov_naive = mean(covered_naive))

# saveRDS(delta_all_ests, file = "delta_all_ests.RDS")
# saveRDS(delta_all_ests_summarized, file = "delta_all_ests_summarized.RDS")


naive_lines <- delta_all_ests_summarized %>%
  distinct(dose_selected, emp_cov_naive) %>%
  mutate(linetype = "Naive Estimator")

cov_plot <- ggplot(delta_all_ests_summarized, aes(x = w, y = emp_cov, color = dose_selected, group = dose_selected)) +
  geom_line(aes(linetype = "Debiased Estimator"), size = 1.5) +
  geom_hline(data = naive_lines,
             aes(yintercept = emp_cov_naive, color = dose_selected, linetype = linetype),
             size = 1.5) +
  scale_linetype_manual(name = "Estimator Type", values = c("Debiased Estimator" = "solid", "Naive Estimator" = "dashed")) +
  labs(color = "Dose Selected", x = "w", y = "Empirical Coverage") +
  ggtitle("Coverage Results from Double Bootstrap Method") +
  geom_hline(yintercept = 0.95, lty = "twodash", col = "blue", lwd = 1.5) +
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "lines"),
    legend.position = "right",
    legend.key.width = unit(1.6, "cm"),  # wider legend lines
    legend.key.height = unit(0.6, "cm"),  # taller legend keys
  )

cov_plot


delta_all_ests %>% filter(w == 1, dose_selected == 2) %>% select(delta_hat) %>% as.vector() %>% unlist() %>%
  density() %>% plot()


lines(density(rnorm(100000, mean = 0.3, sd = 0.1)), , col = "red")
