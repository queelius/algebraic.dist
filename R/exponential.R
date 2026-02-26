#' Construct exponential distribution object.
#'
#' @param rate failure rate
#' @export
exponential <- function(rate) {
  if (!is.numeric(rate) || length(rate) != 1 || is.na(rate) || rate <= 0)
    stop("'rate' must be a positive numeric scalar, got: ", deparse(rate))
  structure(list(rate = rate),
            class = c("exponential", "univariate_dist",
                      "continuous_dist", "dist"))
}

#' Method for obtaining the parameters of an `exponential` object.
#' 
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @export
params.exponential <- function(x) {
  c("rate" = x$rate)
}

#' Function to determine whether an object `x` is an `exponential` object.
#' @param x The object to test
#' @export
is_exponential <- function(x) {
  inherits(x, "exponential")
}

#' Method to obtain the hazard function of an `exponential` object.
#'
#' @param x The `exponential` object to obtain the hazard function of
#' @param ... Additional arguments (not used)
#' @return A function that computes the hazard function of the
#'         exponential distribution at a given point `t` and rate `rate`.
#'         By default, `rate` is the failure rate of object `x`
#'         Also accepts a `log` argument that determines whether
#'         to compute the log of the hazard function.
#' @export
hazard.exponential <- function(x, ...) {
  function(t, rate = x$rate, log = FALSE) {
    stopifnot(is.numeric(rate), rate > 0)

    if (log) {
      ifelse(t <= 0, -Inf, log(rate))
    } else {
      ifelse(t <= 0, 0, rate)
    }
  }
}

#' Method to obtain the pdf of an `exponential` object.
#'
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments (not used)
#' @return A function that computes the pdf of the exponential distribution
#'         at a given point `t`. Also accepts a `rate` argument that
#'         determines the failure rate of the exponential distribution (defaults
#'         to the failure rate of object `x`) and a `log` argument that determines
#'         whether to compute the log of the pdf.
#' @importFrom stats dexp density
#' @export
density.exponential <- function(x, ...) {
  function(t, rate = x$rate, log = FALSE, ...) {
    stopifnot(is.numeric(rate), rate > 0)
    dexp(x = t, rate = rate, log = log, ...)
  }
}

#' Method to sample from an `exponential` object.
#'
#' @param x The `exponential` object to sample from.
#' @param ... Additional arguments to pass (not used)
#' @return A function that allows sampling from the exponential
#'         distribution. Accepts an argument `n` denoting sample
#'         size and `rate` denoting the failure rate. Defaults to
#'         the failure rate of object `x`.
#' @importFrom stats rexp
#' @export
sampler.exponential <- function(x, ...) {
  function(n = 1, rate = x$rate, ...) {
    stopifnot(is.numeric(rate), rate > 0)
    rexp(n = n, rate = rate)
  }
}

#' Method to obtain the mean of an `exponential` object.
#' @param x The `exponential` object to obtain the mean of
#' @param ... Additional arguments (not used)
#' @export
mean.exponential <- function(x, ...) {
  1 / x$rate
}

#' Method to obtain the inverse cdf of an `exponential` object.
#'
#' @param x The object to obtain the inverse cdf of
#' @param ... Additional arguments (not used)
#' @return A function that computes the inverse cdf of the exponential
#'         distribution. Accepts as input a vector `p` probabilities
#'         to compute the inverse cdf, a `rate` value denoting the
#'         failure rate of the exponential distribution, and a logical
#'         `log.p` indicating whether input `p` denotes probability
#'         or log-probability. By default, `rate` is the failure rate
#'         of object `x`.
#' @importFrom stats qexp
#' @export
inv_cdf.exponential <- function(x, ...) {
  function(p, rate = x$rate, lower.tail = TRUE, log.p = FALSE, ...) {
      stopifnot(is.numeric(rate), rate > 0)
      qexp(p = p, rate = rate, lower.tail = lower.tail, log.p = log.p, ...)
  }
}

#' Method to obtain the cdf of an `exponential` object.
#'
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments (not used)
#' @return A function that computes the cdf of the exponential. Accepts as
#'         input a vector `t` at which to compute the cdf, an input `rate`
#'         denoting the failure rate of the exponential distribution, and a
#'         logical `log` indicating whether to compute the log of the cdf.
#'         By default, `rate` is the failure rate of object `x`.
#' @importFrom stats pexp
#' @export
cdf.exponential <- function(x, ...) {
  function(t, rate = x$rate, log.p = FALSE, ...) {
    stopifnot(rate > 0)
    pexp(q = t, rate = rate, log.p = log.p, ...)
  }
}

#' Retrieve the variance of a `exponential` object.
#'
#' @param object The `exponential` object to retrieve the variance for
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `normal` object
#' @export
vcov.exponential <- function(object, ...) {
    1 / object$rate^2
}

#' Support for exponential distribution, the positive real numbers, (0, Inf).
#' @param x The object to obtain the support of
#' @return An `interval` object representing the support of the exponential
#' @export
sup.exponential <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}

#' Method to obtain the cdf of an `exponential` object.
#'
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments (not used)
#' @return A function that computes the cdf of the exponential. Accepts as
#'         input a vector `t` at which to compute the cdf, an input `rate`
#'         denoting the failure rate of the exponential distribution, and a
#'         logical `log` indicating whether to compute the log of the cdf.
#'         By default, `rate` is the failure rate of object `x`.
#' @importFrom stats pexp
#' @export
surv.exponential <- function(x, ...) {
  function(t, rate = x$rate, log.p = FALSE, ...) {
    stopifnot(is.numeric(rate), rate > 0)
    pexp(q = t, rate = rate, log.p = log.p, lower.tail = FALSE, ...)
  }
}

#' Format method for `exponential` objects.
#' @param x The `exponential` object to format
#' @param ... Additional arguments (not used)
#' @return A character string
#' @export
format.exponential <- function(x, ...) {
  sprintf("Exponential distribution (rate = %g)", x$rate)
}

#' Print method for `exponential` objects.
#' @param x The `exponential` object to print.
#' @param ... Additional arguments (not used)
#' @export
print.exponential <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Method to obtain the dimension of an `exponential` object.
#' @param x The `exponential` object to obtain the dimension of
#' @return The dimension of the `exponential` object
#' @export
dim.exponential <- function(x) {
  1
}