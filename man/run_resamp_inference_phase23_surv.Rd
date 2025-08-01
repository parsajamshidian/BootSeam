% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resamp_phase23_surv.R
\name{run_resamp_inference_phase23_surv}
\alias{run_resamp_inference_phase23_surv}
\title{Resampling inference for Survival/Binary Data for seamless Phase 2/3 Trials}
\usage{
run_resamp_inference_phase23_surv(
  n_per_group = c(control = 100, dose1 = 100, dose2 = 100),
  lambda_list = c(control = 2, dose1 = 1, dose2 = 0.5),
  rho = 0.9,
  p_list = c(control = 0.6, dose1 = 0.3, dose2 = 0.15),
  alternative = "less",
  alpha = 0.15,
  alpha_LCL = 0.025,
  ws = seq(0, 1, length.out = 100),
  B = 100
)
}
\arguments{
\item{B}{}
}
\description{
Resampling inference for Survival/Binary Data for seamless Phase 2/3 Trials
}
