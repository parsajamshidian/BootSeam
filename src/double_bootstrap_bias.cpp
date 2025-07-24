#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List double_bootstrap_bias_cpp(NumericVector y0_s1,
                        NumericVector y1_s1,
                        NumericVector y2_s1,
                        int d,
                        int n1,
                        int B1,
                        int B2,
                        Function dose_select,
                        std::string alternative,
                        double alpha,
                        std::string htest_method) {

  NumericVector bias_vals;
  NumericVector delta_hat_nv_star_bs;

  int b1 = 0;

  while (b1 <= B1) {
    IntegerVector boot_ind0 = sample(n1, n1, true) - 1;
    IntegerVector boot_ind1 = sample(n1, n1, true) - 1;
    IntegerVector boot_ind2 = sample(n1, n1, true) - 1;

    NumericVector y0_boot(n1), y1_boot(n1), y2_boot(n1);

    for (int i = 0; i < n1; ++i) {
      y0_boot[i] = y0_s1[boot_ind0[i]];
      y1_boot[i] = y1_s1[boot_ind1[i]];
      y2_boot[i] = y2_s1[boot_ind2[i]];
    }

    List d_sel = dose_select(y0_boot, y1_boot, y2_boot, alternative, alpha, htest_method);
    int d_boot = d_sel["d"];

    if (d_boot != d)
      continue;

    double mean_y0 = mean(y0_boot);
    double mean_y1 = mean(y1_boot);
    double mean_y2 = mean(y2_boot);
    double delta_hat_nv_star_s1 = (d_boot == 1) ? (mean_y1 - mean_y0) : (mean_y2 - mean_y0);

    // Stage 2 of nested bootstrap
    NumericVector bias_star;

    int b2 = 0;
    while (b2 <= B2) {
      IntegerVector ind = sample(n1, n1, true) - 1;
      NumericVector y0_b2(n1), y1_b2(n1), y2_b2(n1);
      for (int i = 0; i < n1; ++i) {
        y0_b2[i] = y0_boot[ind[i]];
        y1_b2[i] = y1_boot[ind[i]];
        y2_b2[i] = y2_boot[ind[i]];
      }

      List d_sel2 = dose_select(y0_b2, y1_b2, y2_b2, alternative, alpha, htest_method);
      int d_boot2 = d_sel2["d"];
      if (d_boot2 != d)
        continue;

      double m_y0 = mean(y0_b2);
      double m_y1 = mean(y1_b2);
      double m_y2 = mean(y2_b2);
      double delta_hat_star_star = (d_boot2 == 1) ? (m_y1 - m_y0) : (m_y2 - m_y0);

      bias_star.push_back(delta_hat_star_star - delta_hat_nv_star_s1);
      ++b2;
    }

    double bias_hat_star = mean(bias_star);
    bias_vals.push_back(delta_hat_nv_star_s1 - (mean_y1 - mean_y0)); // delta_hat_naive_s1 will be passed in R
    delta_hat_nv_star_bs.push_back(delta_hat_nv_star_s1 - bias_hat_star);

    ++b1;
  }

  return List::create(
    Named("bias_vals") = bias_vals,
    Named("delta_hat_nv_star_bs") = delta_hat_nv_star_bs
  );
}
