#' Construct exponential distribution object.
#'
#' @param rate failure rate
#' @export
exp_dist <- function(rate) {
  structure(list(rate = rate),
            class = unique(c("exp_dist",
                             "univariate_dist",
                             "dist")))
}

#' Method for obtaining the parameters of an `exp_dist` object.
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @export
params.exp_dist <- function(x) {
  c("rate" = x$rate)
}

#' Function to determine whether an object `x` is an `exp_dist` object.
#' @export
is_exp_dist <- function(x) {
  inherits(x, "exp_dist")
}

#' Method for obtaining the variance of a `exp_dist` object.
#'
#' @param object The `exp_dist` object to obtain the variance of
#' @param par The rate parameter of the exponential distribution.
#'            If NULL, then the rate of the `exp_dist` object
#'            `object` is used. We view this as the ability to
#'            override its rate, for instance, when we are
#'            performing estimation and the true rate is unknown.
#' @importFrom stats vcov
#' @export
vcov.exp_dist <- function(object, rate = NULL) {
  if (is.null(rate)) {
    rate <- object$rate
  }
  1 / rate^2
}

#' Method to obtain the hazard function of an `exp_dist` object.
#'
#' @param x The `exp_dist` object to obtain the hazard function of
#' @return A function that computes the hazard function of the
#'         exponential distribution at a given point `t` and rate `rate`.
#'         By default, `rate` is the failure rate of object `x`
#'         Also accepts a `log` argument that determines whether
#'         to compute the log of the hazard function.
#' @export
hazard.exp_dist <- function(x) {
  function(t, rate = x$rate, log = FALSE) {
    stopifnot(rate > 0)

    if (log) {
      ifelse(t <= 0, NA, rate)
    } else {
      ifelse(t <= 0, NA, log(rate))
    }
  }
}

#' Method to obtain the pdf of an `exp_dist` object.
#'
#' @param x The object to obtain the pdf of
#' @param rate The rate of the exponential distribution.
#' @return A function that computes the pdf of the exponential distribution
#'         at a given point `t`. Accepts a `log` argument that determines
#'         whether to compute the log of the pdf.
#' 
#' @export
pdf.exp_dist <- function(x) {
  function(t, rate = x$rate, log = FALSE) {
    stopifnot(rate > 0)

    if (log) {
      ifelse(t <= 0, -Inf, log(rate) - rate * t)
    } else {
      ifelse(t <= 0, 0, rate * exp(-rate * t))
    }
  }
}

#' Method to sample from an `exp_dist` object.
#'
#' @param x The `exp_dist` object to sample from.
#' @return A function that allows sampling from the exponential
#'         distribution. Accepts an argument `n` denoting sample
#'         size and `par` denoting the failure rate. If `par` is
#'         NULL, defaults to the failure rate of `x`.
#' @export
sampler.exp_dist <- function(x) {
  function(n = 1, rate = x$rate) {
    stopifnot(rate > 0)
    rexp(n, rate)
  }
}

#' Method to obtain the mean of an `exp_dist` object.
#' @param x The `exp_dist` object to obtain the mean of
#' @param par The rate of the exponential distribution.
#'           If NULL, then the rate of the `exp_dist` object
#'          `x` is used. We view this as the ability to
#'           override its rate, for instance, when we are
#'           performing estimation and the true rate is unknown.
#' @export
mean.exp_dist <- function(x, rate = NULL) {
  if (is.null(rate)) {
    rate <- x$rate
  }
  stopifnot(rate > 0)
  1 / rate
}

#' Method to obtain the inverse cdf of an `exp_dist` object.
#'
#' @param x The object to obtain the inverse cdf of
#' @return A function that computes the inverse cdf of the exponential
#'         distribution. Accepts as input a vector `p` probabilities
#'         to compute the inverse cdf, a `rate` value denoting the
#'         failure rate of the exponential distribution, and a logical
#'         `log.p` indicating whether input `p` denotes probability
#'         or log-probability. By default, `rate` is the failure rate
#'         of object `x`.
#' @export
inv_cdf.exp_dist <- function(x) {
  function(p, rate = x$rate, log.p = FALSE) {
      stopifnot(rate > 0)
      if (log.p) {
        ifelse(p < 0, NA, -log(1-exp(p)) / rate)
      } else {
        ifelse(p < 0 | p > 1, NA, -log(1-p) / rate)
      }
  }
}

#' Method to obtain the cdf of an `exp_dist` object.
#'
#' @param x The object to obtain the pdf of
#' @return A function that computes the cdf of the exponential. Accepts as
#'         input a vector `t` at which to compute the cdf, an input `rate`
#'         denoting the failure rate of the exponential distribution, and a
#'         logical `log` indicating whether to compute the log of the cdf.
#'         By default, `rate` is the failure rate of object `x`.
#' @export
cdf.exp_dist <- function(x) {
  function(t, rate = x$rate, log = FALSE) {
    stopifnot(rate > 0)
    if (log) {
      ifelse(t <= 0, -Inf, log(1 - exp(-rate * t)))
    } else {
      ifelse(t <= 0, 0, 1 - exp(-rate * t))
    }
  }
}

#' Method to obtain the survival function of an `exp_dist` object.
#'
#' @param x The object to obtain the pdf of
#' @return A function that computes the survival function of the exponential,
#'         Accepts as input a vector `t` at which to compute the survival
#'         function, an input `rate` denoting the failure rate of the 
#'         exponential distribution, and a logical `log` indicating whether
#'         to compute the log of the survival function. By default, `rate`
#'         is the failure rate of object `x`.
#' @export
surv.exp_dist <- function(x) {
  function(t, rate = x$rate, log = FALSE) {
    stopifnot(rate > 0)
    if (log) {
      ifelse(t <= 0, 0, -rate * t)
    } else {
      ifelse(t <= 0, 1, exp(-rate * t))
    }
  }
}

#' Support for exponential distribution, the positive real numbers, (0, Inf).
#' @param x The object to obtain the support of
#' @return An `interval` object representing the support of the exponential
#' @export
sup.exp_dist <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}