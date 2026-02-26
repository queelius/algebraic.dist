#' Construct a Poisson distribution object.
#'
#' Creates an S3 object representing a Poisson distribution with rate
#' parameter \eqn{\lambda}.  The PMF is
#' \eqn{P(X = k) = \lambda^k e^{-\lambda} / k!} for \eqn{k = 0, 1, 2, \ldots}.
#'
#' @param lambda Rate parameter (mean), must be a positive scalar.
#' @return A \code{poisson_dist} object with classes
#'   \code{c("poisson_dist", "univariate_dist", "discrete_dist", "dist")}.
#' @export
poisson_dist <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
    stop("'lambda' must be a positive scalar, got: ", deparse(lambda))
  structure(list(lambda = lambda),
            class = c("poisson_dist", "univariate_dist",
                      "discrete_dist", "dist"))
}

#' Test whether an object is a \code{poisson_dist}.
#'
#' @param x The object to test.
#' @return \code{TRUE} if \code{x} inherits from \code{"poisson_dist"},
#'   \code{FALSE} otherwise.
#' @export
is_poisson_dist <- function(x) {
  inherits(x, "poisson_dist")
}

#' Retrieve the parameters of a \code{poisson_dist} object.
#'
#' @param x A \code{poisson_dist} object.
#' @return A named numeric vector with element \code{lambda}.
#' @export
params.poisson_dist <- function(x) {
  c("lambda" = x$lambda)
}

#' Mean of a Poisson distribution.
#'
#' @param x A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return The mean, equal to \code{lambda}.
#' @export
mean.poisson_dist <- function(x, ...) {
  x$lambda
}

#' Variance of a Poisson distribution.
#'
#' For the Poisson distribution, the variance equals \code{lambda}.
#'
#' @param object A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return The variance (scalar), equal to \code{lambda}.
#' @export
vcov.poisson_dist <- function(object, ...) {
  object$lambda
}

#' Dimension of a Poisson distribution (always 1).
#'
#' @param x A \code{poisson_dist} object.
#' @return \code{1}.
#' @export
dim.poisson_dist <- function(x) {
  1
}

#' Format a \code{poisson_dist} object as a character string.
#'
#' @param x A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return A character string describing the distribution.
#' @export
format.poisson_dist <- function(x, ...) {
  sprintf("Poisson distribution (lambda = %g)", x$lambda)
}

#' Print a \code{poisson_dist} object.
#'
#' @param x A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return \code{x}, invisibly.
#' @export
print.poisson_dist <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Support of a Poisson distribution.
#'
#' The Poisson distribution is supported on the non-negative integers
#' \eqn{\{0, 1, 2, \ldots\}}.
#'
#' @param x A \code{poisson_dist} object.
#' @return A \code{countable_set} object with lower bound 0.
#' @export
sup.poisson_dist <- function(x) {
  countable_set$new(lower = 0L)
}

#' Sampler for a Poisson distribution.
#'
#' Returns a function that draws \code{n} independent samples from the
#' Poisson distribution.
#'
#' @param x A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(n = 1, ...)} returning an integer vector
#'   of length \code{n}.
#' @importFrom stats rpois
#' @export
sampler.poisson_dist <- function(x, ...) {
  function(n = 1, ...) {
    rpois(n, lambda = x$lambda)
  }
}

#' Probability mass function for a Poisson distribution.
#'
#' Returns a function that evaluates the Poisson PMF at given points.
#'
#' @param x A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(k, log = FALSE, ...)} returning the
#'   probability mass (or log-probability) at \code{k}.
#' @importFrom stats dpois density
#' @export
density.poisson_dist <- function(x, ...) {
  function(k, log = FALSE, ...) {
    dpois(k, lambda = x$lambda, log = log)
  }
}

#' Cumulative distribution function for a Poisson distribution.
#'
#' Returns a function that evaluates the Poisson CDF at given points.
#'
#' @param x A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(q, log.p = FALSE, ...)} returning the
#'   CDF (or log-CDF) at \code{q}.
#' @importFrom stats ppois
#' @export
cdf.poisson_dist <- function(x, ...) {
  function(q, log.p = FALSE, ...) {
    ppois(q, lambda = x$lambda, log.p = log.p)
  }
}

#' Inverse CDF (quantile function) for a Poisson distribution.
#'
#' Returns a function that computes quantiles of the Poisson distribution.
#'
#' @param x A \code{poisson_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(p, lower.tail = TRUE, log.p = FALSE, ...)}
#'   returning the quantile at probability \code{p}.
#' @importFrom stats qpois
#' @export
inv_cdf.poisson_dist <- function(x, ...) {
  function(p, lower.tail = TRUE, log.p = FALSE, ...) {
    qpois(p, lambda = x$lambda, lower.tail = lower.tail, log.p = log.p)
  }
}

#' Exact expectation for a Poisson distribution.
#'
#' Computes \eqn{E[g(X)]} using truncated summation over the support.
#' The summation is truncated at the \eqn{1 - 10^{-12}} quantile to
#' ensure negligible truncation error.
#'
#' @param x A \code{poisson_dist} object.
#' @param g A function to take the expectation of.
#' @param ... Additional arguments passed to \code{g}.
#' @return The expected value \eqn{E[g(X)]}.
#' @importFrom stats qpois dpois
#' @export
expectation.poisson_dist <- function(x, g, ...) {
  k_max <- qpois(1 - 1e-12, lambda = x$lambda)
  k <- 0:k_max
  sum(g(k, ...) * dpois(k, lambda = x$lambda))
}
