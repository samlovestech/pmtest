#' Corrected Z-test of Looney and Jones
#' @description Perform two sample z test of Looney and Jones on vectors of data.
#' @usage lj.z.test(x, y, alternative = c("two.sided", "less", "greater"))
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
#' lj.z.test(case, control, alternative = "two.sided")
#'
lj.z.test <- function(x, y, alternative = "two.sided") {

  # two samples we are testing shouldn't be empty or with any nonnumeric element
  if(!((all(is.numeric(x)|is.na(x)) & !all(is.na(x))) &
       (all(is.numeric(y)|is.na(y)) & !all(is.na(y))))) {
    stop("Each vector must have at least one numeric values")
  }

  # get paired sample, assign it to t.paired and n.paired
  t.paired <- x[(!is.na(x)) & (!is.na(y))]
  n.paired <- y[(!is.na(x)) & (!is.na(y))]
  n1 <- length(t.paired)  # number of paired sample

  # t.rest is the independent sample from x without NAs
  t.rest <- x[(is.na(y)) & (!is.na(x))]
  n2 <- length(t.rest)

  # n.rest is the independent sample from y without NAs
  n.rest <- y[(is.na(x)) & (!is.na(y))]
  n3 <- length(n.rest)

  #even if n1 == 0, we can't stop the program according the paper.
  if(n1 <= 1) { # if paired sample size is smaller than 1, set the paired sample covariance equals to 0
    paired.cov <- 0
  } else {
    paired.cov <- cov(t.paired,n.paired)
  }

  #calculate mean.star and var.star
  t.mu.star <- mean(c(t.paired,t.rest))
  n.mu.star <- mean(c(n.paired,n.rest))
  t.var.star <- var(c(t.paired,t.rest))
  n.var.star <- var(c(n.paired,n.rest))

  sd <- sqrt(t.var.star/(n1+n2) + n.var.star/(n1+n3) - 2*n1*paired.cov/((n1+n2)*(n1+n3)))
  z.corr = (t.mu.star - n.mu.star)/sd

  if(alternative == "two.sided") {
    p.value <- 2*(1- pnorm(abs(z.corr)))
  }else if(alternative == "greater") {
    p.value <- 1 - pnorm(z.corr)
  }else if(alternative == "less") {
    p.value <- pnorm(z.corr)
  }else{
    stop("arg should be one of \"two.sided\", \"greater\", \"less\"")
  }

  cat("    Corrected Z-test of Looney and Jones\n          p-value is:", p.value, "\n           statistic:", z.corr,"\n")
}
