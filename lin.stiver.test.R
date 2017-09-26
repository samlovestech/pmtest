# MLE based test of Lin and Stivers under heteroscedasticity

#' MLE based test of Lin and Stivers under heteroscedasticity
#' @description Perform two sample MLE based t test of Lin and Stivers under heteroscedasticity on vectors of data.
#' @usage lin.stivers.test(x, y, alternative = c("two.sided", "less", "greater"))
#' @param x a (non-empty) numeric vector of data values.
#' @param x a (non-empty) numeric vector of data values.
#' @param alternative	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @details
#' The formula interface is only applicable for the 2-sample tests.
#' alternative = "greater" is the alternative that `x` has a larger mean than `y`.`x` and `y` must be the same length.
#' Missing values are treated such that matched pairs are grouped, and non matches in `x` and `y` are treated as independent samples.
#' @return A list containing the following components:
#' \item{statistic}{the value of t-statistic}
#' \item{p.value}{p-value of the test}
#' \item{df}{degree of freedom}
#'
#' @examples
#' case <- c(10,NA,12,14,2,4)
#' control <- c(8,4,NA,13,4,5)
#' lin.stivers.test(case, control, alternative = "two.sided")
#'
lin.stivers.test <- function(x, y, alternative = "two.sided") {

  # two samples we are testing shouldn't be empty or with any nonnumeric element
  if(!((all(is.numeric(x)|is.na(x)) & !all(is.na(x))) &
       (all(is.numeric(y)|is.na(y)) & !all(is.na(y))))) {
    stop("Each vector must have at least one numeric values")
  }

  # get paired sample, assign it to t.paired and n.paired
  t.paired <- x[(!is.na(x)) & (!is.na(y))]
  n.paired <- y[(!is.na(x)) & (!is.na(y))]
  n1 <- length(t.paired)  # number of paired sample
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

  # if the variance is equal at level 0.05, we give warning.
  x.na.rm <- x[!is.na(x)]
  y.na.rm <- y[!is.na(y)]
  f.test = var.test(x.na.rm, y.na.rm)

  r <- var(t.paired, n.paired)/(sd(t.paired)*sd(n.paired))
  f <- n1*(n1+n3+n2*cov(t.paired, n.paired)/var(t.paired))/((n1+n2)*(n1+n3) - n2*n3*r*r)
  g <- n1*(n1+n2+n3*cov(t.paired, n.paired)/var(n.paired))/((n1+n2)*(n1+n3) - n2*n3*r*r)

  if(n2 == 0 & n3 != 0) {
    t.rest <- 0
    v1.numerator <- (f^2/n1)*var(t.paired)*(n1-1) + (g^2/n1 + (1-g)^2/n3)*var(n.paired)*(n1-1) - 2*f*g*cov(t.paired,n.paired)*(n1-1)/n1
  } else if (n3 == 0 & n2 != 0) {
    n.rest <- 0
    v1.numerator <- (f^2/n1 + (1-f)^2/n2)*var(t.paired)*(n1-1) + (g^2/n1)*var(n.paired)*(n1-1) - 2*f*g*cov(t.paired,n.paired)*(n1-1)/n1
  } else if(n2 == 0 & n3 == 0) {
    t.rest <- 0
    n.rest <- 0
    v1.numerator <- (f^2/n1)*var(t.paired)*(n1-1) + (g^2/n1)*var(n.paired)*(n1-1) - 2*f*g*cov(t.paired,n.paired)*(n1-1)/n1
  } else {
    v1.numerator <- (f^2/n1 + (1-f)^2/n2)*var(t.paired)*(n1-1) + (g^2/n1 + (1-g)^2/n3)*var(n.paired)*(n1-1) - 2*f*g*cov(t.paired,n.paired)*(n1-1)/n1
  }

  v1 <- v1.numerator/(n1-1)
  diff <- f*(mean(t.paired) - mean(t.rest)) - g*(mean(n.paired) - mean(n.rest)) + mean(t.rest) - mean(n.rest)
  z.ls <- diff/sqrt(v1)
  df <- n1

  if(alternative == "two.sided") {
    p.value <- 2*(1- pt(abs(z.ls), df))
  }else if(alternative == "greater") {
    p.value <- 1 - pt(z.ls, df)
  }else if(alternative == "less") {
    p.value <- pt(z.ls, df)
  }else {
    stop("arg should be one of \"two.sided\", \"greater\", \"less\"")
  }

  # relative warnings, should put it after alternative = typo, since if typo happens, i dont want warning happens.
  if(f.test$p.value <= 0.05) {
    warning("Warning: The assumption of lin stivers test is the different variance of two samples. Since the variance are the same between two samples at α = 0.05, ekbohm test is recommended", noBreaks. = T)
  }
  # if r ~ 1, it will cause problems
  if(r >= 0.99 ) {
    warning("the correlation of paried sample is close to 1, this lin stiver test might not be accurate",noBreaks. = T)
  }
  if(n2 <= 1 | n3 <= 1) {
    warning("since there is only few NA in two samples, this lin stiver test might not be accurate.",noBreaks. = T)
  }

  cat("    MLE based test of Lin and Stivers under heteroscedasticity\n               p-value is:", p.value, "\n                statistic:", z.ls, "\n                       df:", df,"\n")
}
