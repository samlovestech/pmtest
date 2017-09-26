# MLE based test of Ekbohm under homoscedasticity
# when the variances of tumor and normal are equal, Ekbohm [5] suggested the following MLE based test statistic

#' MLE based test of Ekbohm under homoscedasticity
#' @description Perform two sample MLE based t test of Ekbohm under homoscedasticity on vectors of data.
#' @usage ekbohm.test(x, y, alternative = c("two.sided", "less", "greater"))
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
#' \item{df}{degree of freedom}
#'
#' @examples
#' case <- c(10,NA,12,14,2,4)
#' control <- c(8,4,NA,13,4,5)
#' ekbohm.test(case, control, alternative = "two.sided")
#'
ekbohm.test <- function(x, y, alternative = "two.sided") {

  # two samples we are testing shouldn't be empty or with any nonnumeric element
  if(!((all(is.numeric(x)|is.na(x)) & !all(is.na(x))) &
       (all(is.numeric(y)|is.na(y)) & !all(is.na(y))))) {
    stop("Each vector must have at least one numeric values")
  }

  # get the paired sample, assign it to t.paired and n.paired
  t.paired <- x[(!is.na(x)) & (!is.na(y))]
  n.paired <- y[(!is.na(x)) & (!is.na(y))]
  n1 <- length(t.paired)  # number of paired sample
  # if paired sample size is smaller than 1, stop this testing.
  if(n1 <= 1) {
    stop("not enough paired observations")
  }

  # get the rest sample which is in T(not NA) and in N(is NA)
  t.rest <- x[(is.na(y)) & (!is.na(x))]
  n2 <- length(t.rest)

  # get the rest sample which is in N(not NA) and in T(is NA)
  n.rest <- y[(is.na(x)) & (!is.na(y))]
  n3 <- length(n.rest)

  # Since this Ekbohm test is based on the same variance, we give warning if two samples have different variance at level 0.05
  x.na.rm <- x[!is.na(x)]
  y.na.rm <- y[!is.na(y)]
  f.test = var.test(x.na.rm, y.na.rm)

  if(n2 <= 1) {# independent sample from x
    if(n2 ==0) {
      t.rest <- 0
    }
    var.t.rest <- 0
  } else{
    var.t.rest <- var(t.rest)
  }

  if(n3 <= 1) {# independent sample from y
    if(n3 == 0) {
      n.rest <- 0
    }
    var.n.rest <- 0
  } else{
    var.n.rest <- var(n.rest)
  }

  r <- cov(t.paired, n.paired)/(sd(t.paired)*sd(n.paired))
  f <- n1*(n1+n3+n2*r)/((n1+n2)*(n1+n3) - n2*n3*r^2)
  g <- n1*(n1+n2+n3*r)/((n1+n2)*(n1+n3) -n2*n3*r^2)
  var.hat <- (var(t.paired)*(n1-1) + var(n.paired)*(n1-1) + (1+r^2)*(var.t.rest*(n2-1) + var.n.rest*(n3-1)))/(2*(n1-1) + (1+r^2)*(n2+n3-2))
  v1.star <- var.hat*((2*n1*(1-r) + (n2+n3)*(1-r^2))/((n1+n2)*(n1+n3) - n2*n3*r^2))

  diff <- f*(mean(t.paired) - mean(t.rest)) - g*(mean(n.paired) - mean(n.rest)) + mean(t.rest) - mean(n.rest)
  z.e <- diff/sqrt(v1.star)
  df <- n1

  if(alternative == "two.sided") {
    p.value <- 2*(1 - pt(abs(z.e), df))
  }else if(alternative == "greater") {
    p.value <- 1 - pt(z.e, df)
  }else if(alternative == "less") {
    p.value <- pt(z.e, df)
  }else {
    stop("arg should be one of \"two.sided\", \"greater\", \"less\"")
  }

  # relative warnings, should put it after alternative = typo, since if typo happens, i dont want warning happens.
  if(f.test$p.value >= 0.05) {
    warning("Warning: The assumption of Ekbohm test is the same variance of two samples. Since the variance are different between two samples at α = 0.05, Lin Stiver test is recommended", noBreaks. = T)
  }
  #  if r ~ 1, it will cause problems
  if(r >= 0.99 ) {
    warning("the correlation of paried sample is close to 1, this ekbohm test might not be accurate", noBreaks. = T)
  }
  if(n2 <= 1 | n3 <= 1) {
    warning("since there is only few NA in two samples, this ekbohm test might not be accurate.",noBreaks. = T)
  }

  cat("    MLE based test of Ekbohm under homoscedasticity\n            p-value is:", p.value, "\n             statistic:", z.e, "\n                    df:", df,"\n")
}
