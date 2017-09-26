#' Modified t-statistic of Kim et al
#' @description Perform two sample Kim t test on vectors of data.
#' @usage kim.t.test(x, y, alternative = c("two.sided", "less", "greater"))
#' @param x a (non-empty) numeric vector of data values.
#' @param y a (non-empty) numeric vector of data values.
#' @param alternative	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @details
#' The formula interface is only applicable for the 2-sample tests.
#' alternative = "greater" is the alternative that `x` has a larger mean than `y`.`x` and `y` must be the same length.
#' Missing values are treated such that matched pairs are grouped, and non matches in `x` and `y` are treated as independent samples.
#'
#' @return A list containing the following components:
#' \item{statistic}{the value of t-statistic}
#' \item{p.value}{p-value of the test}
#'
#' @examples
#' case <- c(10,NA,12,14,2,4)
#' control <- c(8,4,NA,13,4,5)
#' kim.t.test(case, control, alternative = "two.sided")
#'
#'
kim.t.test <- function(x, y, alternative = "two.sided") {

  # two samples we are testing shouldn't be empty or with any nonnumeric element
  if(!((all(is.numeric(x)|is.na(x)) & !all(is.na(x))) &
       (all(is.numeric(y)|is.na(y)) & !all(is.na(y))))) {
    stop("Each vector must have at least one numeric values")
  }

  # x is the first sample, y is the second sample
  # d is the difference between the paired samples where there are no NA values
  d <- x-y
  d <- d[!is.na(d)]
  n1 <- length(d)
  # if paired sample size is smaller than 1, stop this testing.
  if(n1 <= 1) {
    stop("not enough paired observations")
  }

  # t.rest is the independent sample from x without NAs
  t.rest <- x[(is.na(y)) & (!is.na(x))]
  n2 <- length(t.rest)

  # n.rest is the independent sample from y without NAs
  n.rest <- y[(is.na(x)) & (!is.na(y))]
  n3 <- length(n.rest)


  if(n2 <= 1) { # independent sample from x
    if(n2 ==0) {
      t.rest <- 0
    }
    var.t.rest <- 0
  } else {
    var.t.rest <- var(t.rest)
  }
  if(n3 <= 1) {# independent sample from y
    if(n3 == 0) {
      n.rest <- 0
    }
    var.n.rest <- 0
  } else {
    var.n.rest <- var(n.rest)
  }


  # as the paper says the null distribution of t is approximated with N(0,1), so we use z-table
  # calculate statistic t3
  n.h <- 2/(1/n2 + 1/n3)  # get n(h): harmonic mean of t.num and n.num
  if(n3 == 0 | n2 == 0) {
    sd <- sqrt(n1*var(d))
    t3 <- n1*mean(d)/sd
  } else {
    sd <- sqrt(n1*var(d) + n.h^2*(var.n.rest/n3 + var.t.rest/n2))
    t3 <- (n1*mean(d) + n.h*(mean(t.rest) - mean(n.rest)))/sd
  }

  # p.value
  if(alternative == "two.sided") {
    p.value <- 2*(1- pnorm(abs(t3)))
  }else if(alternative == "greater") {
    p.value <- 1 - pnorm(t3)
  }else if(alternative == "less") {
    p.value <- pnorm(t3)
  }else {
    stop("arg should be one of \"two.sided\", \"greater\", \"less\"")
  }

  # if(n2 <= 1 | n3 <= 1) {
  #   warning("since there is only few NA in two samples, this kim.t.test might not be accurate.",noBreaks. = T)
  # }

  cat("    Modified t test of Kim et al\n      p-value is:", p.value, "\n       statistic:", t3,"\n")
}

