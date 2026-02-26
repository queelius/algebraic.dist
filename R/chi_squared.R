#' Construct a chi-squared distribution object.
#'
#' @param df Degrees of freedom (positive scalar)
#' @return A `chi_squared` object
#' @export
chi_squared <- function(df) {
  if (!is.numeric(df) || length(df) != 1 || is.na(df) || df <= 0)
    stop("'df' must be a positive scalar, got: ", deparse(df))
  structure(list(df = df),
            class = c("chi_squared", "univariate_dist",
                      "continuous_dist", "dist"))
}

#' Test whether an object is a `chi_squared`.
#' @param x The object to test
#' @return Logical; TRUE if `x` inherits from `chi_squared`
#' @export
is_chi_squared <- function(x) {
  inherits(x, "chi_squared")
}

#' Method for obtaining the parameters of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @return A named numeric vector of parameters
#' @export
params.chi_squared <- function(x) {
  c(df = x$df)
}

#' Retrieve the mean of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return The mean, equal to `df`
#' @export
mean.chi_squared <- function(x, ...) {
  x$df
}

#' Retrieve the variance of a `chi_squared` object.
#'
#' @param object The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return The variance, `2 * df`
#' @export
vcov.chi_squared <- function(object, ...) {
  2 * object$df
}

#' Retrieve the dimension of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @return 1 (univariate)
#' @export
dim.chi_squared <- function(x) {
  1
}

#' Format a `chi_squared` object as a character string.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A character string describing the distribution
#' @export
format.chi_squared <- function(x, ...) {
  sprintf("Chi-squared distribution (df = %g)", x$df)
}

#' Print method for `chi_squared` objects.
#'
#' @param x The `chi_squared` object to print
#' @param ... Additional arguments (not used)
#' @export
print.chi_squared <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Method for sampling from a `chi_squared` object.
#'
#' @param x The `chi_squared` object to sample from
#' @param ... Additional arguments (not used)
#' @return A function that generates `n` samples from the chi-squared
#'         distribution
#' @importFrom stats rchisq
#' @export
sampler.chi_squared <- function(x, ...) {
  function(n = 1) {
    rchisq(n = n, df = x$df)
  }
}

#' Method for obtaining the density (pdf) of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A function that computes the pdf at point(s) `t`
#' @importFrom stats dchisq density
#' @export
density.chi_squared <- function(x, ...) {
  function(t, log = FALSE) {
    dchisq(x = t, df = x$df, log = log)
  }
}

#' Method for obtaining the cdf of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A function that computes the cdf at point(s) `t`
#' @importFrom stats pchisq
#' @export
cdf.chi_squared <- function(x, ...) {
  function(t, log.p = FALSE) {
    pchisq(q = t, df = x$df, log.p = log.p)
  }
}

#' Method for obtaining the inverse cdf (quantile function) of a `chi_squared`
#' object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A function that computes the quantile at probability `p`
#' @importFrom stats qchisq
#' @export
inv_cdf.chi_squared <- function(x, ...) {
  function(p, lower.tail = TRUE, log.p = FALSE) {
    qchisq(p = p, df = x$df, lower.tail = lower.tail, log.p = log.p)
  }
}

#' Support for chi-squared distribution, the positive real numbers (0, Inf).
#'
#' @param x The `chi_squared` object
#' @return An `interval` object representing (0, Inf)
#' @export
sup.chi_squared <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}

#' Method for obtaining the survival function of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A function that computes S(t) = P(X > t)
#' @export
surv.chi_squared <- function(x, ...) {
  function(t, log.p = FALSE) {
    pchisq(q = t, df = x$df, lower.tail = FALSE, log.p = log.p)
  }
}

#' Method for obtaining the hazard function of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A function that computes h(t) = f(t) / S(t)
#' @export
hazard.chi_squared <- function(x, ...) {
  f <- density(x)
  S <- surv(x)
  function(t, log = FALSE) {
    h <- f(t) / S(t)
    if (log) log(h) else h
  }
}
