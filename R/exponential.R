#' Construct exponential distribution object.
#'
#' @param rate failure rate
#' @return An \code{exponential} distribution object.
#' @examples
#' x <- exponential(rate = 2)
#' mean(x)
#' vcov(x)
#' format(x)
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
#' @examples
#' x <- exponential(rate = 0.5)
#' params(x)
#' @export
params.exponential <- function(x) {
  c("rate" = x$rate)
}

#' Function to determine whether an object `x` is an `exponential` object.
#' @param x The object to test
#' @return Logical; \code{TRUE} if \code{x} is an \code{exponential}.
#' @examples
#' is_exponential(exponential(1))
#' is_exponential(normal(0, 1))
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
#' @examples
#' x <- exponential(rate = 2)
#' h <- hazard(x)
#' h(1)
#' h(5)
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
#' @examples
#' x <- exponential(rate = 2)
#' f <- density(x)
#' f(0)
#' f(1)
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
#' @examples
#' x <- exponential(rate = 2)
#' s <- sampler(x)
#' set.seed(42)
#' s(5)
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
#' @return The mean of the exponential distribution (\code{1 / rate}).
#' @examples
#' x <- exponential(rate = 0.5)
#' mean(x)
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
#' @examples
#' x <- exponential(rate = 1)
#' q <- inv_cdf(x)
#' q(0.5)
#' q(0.95)
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
#' @examples
#' x <- exponential(rate = 1)
#' F <- cdf(x)
#' F(1)
#' F(2)
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
#' @examples
#' x <- exponential(rate = 2)
#' vcov(x)
#' @export
vcov.exponential <- function(object, ...) {
    1 / object$rate^2
}

#' Support for exponential distribution, the positive real numbers, (0, Inf).
#' @param x The object to obtain the support of
#' @return An `interval` object representing the support of the exponential
#' @examples
#' x <- exponential(rate = 1)
#' sup(x)
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
#' @examples
#' x <- exponential(rate = 1)
#' S <- surv(x)
#' S(1)
#' S(2)
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
#' @examples
#' format(exponential(rate = 2))
#' @export
format.exponential <- function(x, ...) {
  sprintf("Exponential distribution (rate = %g)", x$rate)
}

#' Print method for `exponential` objects.
#' @param x The `exponential` object to print.
#' @param ... Additional arguments (not used)
#' @return \code{x}, invisibly.
#' @examples
#' print(exponential(rate = 2))
#' @export
print.exponential <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Method to obtain the dimension of an `exponential` object.
#' @param x The `exponential` object to obtain the dimension of
#' @return The dimension of the `exponential` object
#' @examples
#' dim(exponential(rate = 1))
#' @export
dim.exponential <- function(x) {
  1
}