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
#' 
#' If `par` is NULL, then the rate of the `exp_dist` object `object` is used.
#' We view this as the ability to override its rate, for instance, when we are
#' performing estimation and the true rate is unknown.
#'
#' @param x The object to obtain the parameters of
#' @param par The rate parameter of the exponential distribution.
#' @return A named vector of parameters
#' @export
params.exp_dist <- function(x, par = NULL) {
  if (is.null(par)) {
    par <- x$par
  }
  c("rate" = par)
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
#' @importFrom stats vcov
#' @export
vcov.exp_dist <- function(object, par = NULL) {
  par <- params(object, par)
  1 / par^2
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
  function(t, par = x$rate, log = FALSE) {
    stopifnot(par > 0)

    if (log) {
      ifelse(t <= 0, -Inf, log(par))
    } else {
      ifelse(t <= 0, 0, par)
    }
  }
}

#' Method to obtain the pdf of an `exp_dist` object.
#'
#' @param x The object to obtain the pdf of
#' @param par The rate of the exponential distribution.
#' @return A function that computes the pdf of the exponential distribution
#'         at a given point `t`. Accepts a `log` argument that determines
#'         whether to compute the log of the pdf.
#' 
#' @export
pdf.exp_dist <- function(x) {
  function(t, par = x$rate, log = FALSE) {
    stopifnot(par > 0)
    dexp(t, par, log)
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
  function(n = 1, par = x$rate) {
    stopifnot(par > 0)
    rexp(n, par)
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
mean.exp_dist <- function(x, par = NULL) {
  if (is.null(par)) {
    par <- x$rate
  }
  stopifnot(par > 0)
  1 / par
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
inv_cdf.exp_dist <- function(x, rate = NULL, ...) {
  function(p, rate = x$rate, lower.tail = TRUE, log.p = FALSE) {
      stopifnot(rate > 0)
      qexp(p = p, rate = rate, lower.tail = lower.tail, log.p = log.p)
  }
}

#' Method to obtain the cdf of an `exp_dist` object.
#'
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments (not used)
#' @return A function that computes the cdf of the exponential. Accepts as
#'         input a vector `t` at which to compute the cdf, an input `rate`
#'         denoting the failure rate of the exponential distribution, and a
#'         logical `log` indicating whether to compute the log of the cdf.
#'         By default, `rate` is the failure rate of object `x`.
#' @export
cdf.exp_dist <- function(x, ...) {
  function(t, rate = x$rate, lower.tail = TRUE, log.p = FALSE) {
    stopifnot(rate > 0)
    pexp(q = t, rate = rate, lower.tail = lower.tail, log.p = log.p)
  }
}

#' Support for exponential distribution, the positive real numbers, (0, Inf).
#' @param x The object to obtain the support of
#' @return An `interval` object representing the support of the exponential
#' @export
sup.exp_dist <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}

#' Method to obtain the cdf of an `exp_dist` object.
#'
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments (not used)
#' @return A function that computes the cdf of the exponential. Accepts as
#'         input a vector `t` at which to compute the cdf, an input `rate`
#'         denoting the failure rate of the exponential distribution, and a
#'         logical `log` indicating whether to compute the log of the cdf.
#'         By default, `rate` is the failure rate of object `x`.
#' @export
surv.exp_dist <- function(x, ...) {
  function(t, rate = x$rate, log.p = FALSE) {
    stopifnot(rate > 0)
    pexp(q = t, rate = rate, log.p = log.p)
  }
}

#' This is the expectation of a function `g` with respect to an
#' exponential distribution `x` of type `exp_dist`.
#'
#' @param x The disrtibution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass.
#' @export
expectation.exp_dist <- function(x, g, rate = NULL, ...) {

  rate <- params(x, rate)
  stopifnot(is.function(g), rate > 0)
  f <- pdf(x)
  integrate(function(t) g(t) * f(t, ...), 0, Inf)
}


#' Print method for `exp_dist` objects.
#' @param x The `exp_dist` object to print.
#' @param ... Additional arguments (not used)
#' @export
print.exp_dist <- function(x, ...) {
  cat("Exponential distribution with failure rate", x$rate, "\n")
}