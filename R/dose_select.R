# Dose selection algorithm
#'
#'
#' @param y0
#' @param y1
#' @param y2
#' @param alternative
#' @param alpha
#' @param htest_method
#'
#' @returns
#' @export
#'
#' @examples
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
