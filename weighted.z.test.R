#' Weighted Z-test combination
#' @description Perform two sample weighted z test on vectors of data.
#' @usage weighted.z.test(x, y, alternative = c("two.sided", "less", "greater"))
#' @param x a (non-empty) numeric vector of data values.
#' @param y a (non-empty) numeric vector of data values.
#' @param alternative	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @details
#' The formula interface is only applicable for the 2-sample tests.
#' alternative = "greater" is the alternative that `x` has a larger mean than `y`.`x` and `y` must be the same length.
#' Missing values are treated such that matched pairs are grouped, and non matches in `x` and `y` are treated as independent samples.
#' @return A list containing the following components:
#' \item{statistic}{the value of t-statistic}
#' \item{p.value}{p-value of the test}
#'
#' @examples
#' case <- c(10,NA,12,14,2,4)
#' control <- c(8,4,NA,13,4,5)
#' weighted.z.test(case, control, alternative = "two.sided")
#'
weighted.z.test <- function(x, y, alternative = "greater") {

  # two samples we are testing shouldn't be empty or with any nonnumeric element
  if(!((all(is.numeric(x)|is.na(x)) & !all(is.na(x))) &
       (all(is.numeric(y)|is.na(y)) & !all(is.na(y))))) {
    stop("Each vector must have at least one numeric values")
  }

  # get paired sample, assign it to t.paired and n.paired
  t.paired <- x[(!is.na(x)) & (!is.na(y))]
  n.paired <- y[(!is.na(x)) & (!is.na(y))]
  n1 <- length(t.paired)  # number of paired sample
  if(n1 <= 1) {
    stop("not enough paired observations")
  }

  # t.rest is the independent sample from x without NAs
  t.rest <- x[(is.na(y)) & (!is.na(x))]
  n2 <- length(t.rest)

  # n.rest is the independent sample from y without NAs
  n.rest <- y[(is.na(x)) & (!is.na(y))]
  n3 <- length(n.rest)

  if(n2 <= 1 | n3 <= 1) {
    p.value <- t.test(t.paired, n.paired, paired = T,alternative = alternative)$p.value
    cat("    Weighted Z-test Combination\n      p-value is:",p.value)
    stop("Since there's only few NA value in the paired sample, the result is given by paired t.test ")
  }

  p.1i <- t.test(t.paired, n.paired, paired = T, alternative = "greater")$p.value
  p.2i <- t.test(t.rest, n.rest, alternative = "greater")$p.value
  z.1i <- qnorm(1-p.1i)
  z.2i <- qnorm(1-p.2i)
  w1 <- sqrt(2*n1)
  w2 <- sqrt(n2+n3)  # by Zaykin[10]... square root of the sample sizes in practice
  p.ci <- 1 - pnorm((w1*z.1i + w2*z.2i)/sqrt(w1^2 + w2^2))

  if(alternative == "two.sided") {
    if(p.ci < 0.5) {
      p.value <- 2*p.ci
    }else {
      p.value <- 2*(1 - p.ci)
    }
  }else if(alternative == "greater") {
    p.value <- p.ci
  }else if(alternative == "less") {
    p.value <- 1 - p.ci
  }else {
    stop("arg should be one of \"two.sided\", \"greater\", \"less\"")
  }

  cat("    Weighted Z-test Combination\n      p-value is:",p.value,"\n")
}


