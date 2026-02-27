#' Construct a chi-squared distribution object.
#'
#' @param df Degrees of freedom (positive scalar)
#' @return A `chi_squared` object
#' @examples
#' x <- chi_squared(df = 5)
#' mean(x)
#' vcov(x)
#' format(x)
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
#' @examples
#' is_chi_squared(chi_squared(3))
#' is_chi_squared(normal(0, 1))
#' @export
is_chi_squared <- function(x) {
  inherits(x, "chi_squared")
}

#' Method for obtaining the parameters of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @return A named numeric vector of parameters
#' @examples
#' params(chi_squared(5))
#' @export
params.chi_squared <- function(x) {
  c(df = x$df)
}

#' Retrieve the mean of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return The mean, equal to `df`
#' @examples
#' mean(chi_squared(10))
#' @export
mean.chi_squared <- function(x, ...) {
  x$df
}

#' Retrieve the variance of a `chi_squared` object.
#'
#' @param object The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return The variance, `2 * df`
#' @examples
#' vcov(chi_squared(10))
#' @export
vcov.chi_squared <- function(object, ...) {
  2 * object$df
}

#' Retrieve the dimension of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @return 1 (univariate)
#' @examples
#' dim(chi_squared(5))
#' @export
dim.chi_squared <- function(x) {
  1
}

#' Format a `chi_squared` object as a character string.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A character string describing the distribution
#' @examples
#' format(chi_squared(5))
#' @export
format.chi_squared <- function(x, ...) {
  sprintf("Chi-squared distribution (df = %g)", x$df)
}

#' Print method for `chi_squared` objects.
#'
#' @param x The `chi_squared` object to print
#' @param ... Additional arguments (not used)
#' @return \code{x}, invisibly.
#' @examples
#' print(chi_squared(5))
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
#' @examples
#' x <- chi_squared(5)
#' s <- sampler(x)
#' set.seed(42)
#' s(5)
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
#' @examples
#' x <- chi_squared(5)
#' f <- density(x)
#' f(5)
#' f(10)
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
#' @examples
#' x <- chi_squared(5)
#' F <- cdf(x)
#' F(5)
#' F(10)
#' @importFrom stats pchisq
#' @export
cdf.chi_squared <- function(x, ...) {
  function(q, log.p = FALSE) {
    pchisq(q = q, df = x$df, log.p = log.p)
  }
}

#' Method for obtaining the inverse cdf (quantile function) of a `chi_squared`
#' object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A function that computes the quantile at probability `p`
#' @examples
#' x <- chi_squared(5)
#' q <- inv_cdf(x)
#' q(0.5)
#' q(0.95)
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
#' @examples
#' sup(chi_squared(5))
#' @export
sup.chi_squared <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}

#' Method for obtaining the survival function of a `chi_squared` object.
#'
#' @param x The `chi_squared` object
#' @param ... Additional arguments (not used)
#' @return A function that computes S(t) = P(X > t)
#' @examples
#' x <- chi_squared(5)
#' S <- surv(x)
#' S(5)
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
#' @examples
#' x <- chi_squared(5)
#' h <- hazard(x)
#' h(5)
#' @export
hazard.chi_squared <- function(x, ...) {
  function(t, log = FALSE) {
    log_f <- dchisq(t, df = x$df, log = TRUE)
    log_S <- pchisq(t, df = x$df, lower.tail = FALSE, log.p = TRUE)
    log_h <- log_f - log_S
    if (log) log_h else exp(log_h)
  }
}
