#' Construct a beta distribution object.
#'
#' Creates an S3 object representing a beta distribution with shape
#' parameters \code{shape1} and \code{shape2}.  The PDF on \eqn{(0, 1)} is
#' \deqn{f(x) = \frac{x^{a-1}(1-x)^{b-1}}{B(a,b)}}
#' where \eqn{a} = \code{shape1}, \eqn{b} = \code{shape2}, and \eqn{B(a,b)}
#' is the beta function.
#'
#' @param shape1 First shape parameter, must be a positive scalar.
#' @param shape2 Second shape parameter, must be a positive scalar.
#' @return A \code{beta_dist} object with classes
#'   \code{c("beta_dist", "univariate_dist", "continuous_dist", "dist")}.
#' @examples
#' x <- beta_dist(shape1 = 2, shape2 = 5)
#' mean(x)
#' vcov(x)
#' format(x)
#' @export
beta_dist <- function(shape1, shape2) {
  if (!is.numeric(shape1) || length(shape1) != 1 || is.na(shape1) || shape1 <= 0)
    stop("'shape1' must be a positive scalar, got: ", deparse(shape1))
  if (!is.numeric(shape2) || length(shape2) != 1 || is.na(shape2) || shape2 <= 0)
    stop("'shape2' must be a positive scalar, got: ", deparse(shape2))
  structure(list(shape1 = shape1, shape2 = shape2),
            class = c("beta_dist", "univariate_dist",
                      "continuous_dist", "dist"))
}

#' Test whether an object is a \code{beta_dist}.
#'
#' @param x The object to test.
#' @return \code{TRUE} if \code{x} inherits from \code{"beta_dist"},
#'   \code{FALSE} otherwise.
#' @examples
#' is_beta_dist(beta_dist(2, 5))
#' is_beta_dist(normal(0, 1))
#' @export
is_beta_dist <- function(x) {
  inherits(x, "beta_dist")
}

#' Retrieve the parameters of a \code{beta_dist} object.
#'
#' @param x A \code{beta_dist} object.
#' @return A named numeric vector with elements \code{shape1} and \code{shape2}.
#' @examples
#' params(beta_dist(2, 5))
#' @export
params.beta_dist <- function(x) {
  c(shape1 = x$shape1, shape2 = x$shape2)
}

#' Mean of a beta distribution.
#'
#' Computes \eqn{\alpha / (\alpha + \beta)} where \eqn{\alpha} = \code{shape1}
#' and \eqn{\beta} = \code{shape2}.
#'
#' @param x A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return The mean of the distribution.
#' @examples
#' mean(beta_dist(2, 5))
#' @export
mean.beta_dist <- function(x, ...) {
  x$shape1 / (x$shape1 + x$shape2)
}

#' Variance of a beta distribution.
#'
#' Computes \eqn{\alpha\beta / ((\alpha+\beta)^2 (\alpha+\beta+1))}.
#'
#' @param object A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return The variance (scalar).
#' @examples
#' vcov(beta_dist(2, 5))
#' @export
vcov.beta_dist <- function(object, ...) {
  a <- object$shape1
  b <- object$shape2
  a * b / ((a + b)^2 * (a + b + 1))
}

#' Dimension of a beta distribution (always 1).
#'
#' @param x A \code{beta_dist} object.
#' @return \code{1}.
#' @examples
#' dim(beta_dist(2, 5))
#' @export
dim.beta_dist <- function(x) {
  1
}

#' Format a \code{beta_dist} object as a character string.
#'
#' @param x A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return A character string describing the distribution.
#' @examples
#' format(beta_dist(2, 5))
#' @export
format.beta_dist <- function(x, ...) {
  sprintf("Beta distribution (shape1 = %g, shape2 = %g)", x$shape1, x$shape2)
}

#' Print a \code{beta_dist} object.
#'
#' @param x A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return \code{x}, invisibly.
#' @examples
#' print(beta_dist(2, 5))
#' @export
print.beta_dist <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Support of a beta distribution.
#'
#' The beta distribution is supported on the open interval \eqn{(0, 1)}.
#'
#' @param x A \code{beta_dist} object.
#' @return An \code{interval} object representing \eqn{(0, 1)}.
#' @examples
#' sup(beta_dist(2, 5))
#' @export
sup.beta_dist <- function(x) {
  interval$new(lower = 0, upper = 1,
               lower_closed = FALSE, upper_closed = FALSE)
}

#' Sampler for a beta distribution.
#'
#' Returns a function that draws \code{n} independent samples from the
#' beta distribution.
#'
#' @param x A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(n = 1, ...)} returning a numeric vector
#'   of length \code{n}.
#' @examples
#' x <- beta_dist(2, 5)
#' s <- sampler(x)
#' set.seed(42)
#' s(5)
#' @importFrom stats rbeta
#' @export
sampler.beta_dist <- function(x, ...) {
  function(n = 1, ...) {
    rbeta(n, shape1 = x$shape1, shape2 = x$shape2)
  }
}

#' Probability density function for a beta distribution.
#'
#' Returns a function that evaluates the beta PDF at given points.
#'
#' @param x A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log = FALSE, ...)} returning the
#'   density (or log-density) at \code{t}.
#' @examples
#' x <- beta_dist(2, 5)
#' f <- density(x)
#' f(0.3)
#' f(0.5)
#' @importFrom stats dbeta density
#' @export
density.beta_dist <- function(x, ...) {
  function(t, log = FALSE, ...) {
    dbeta(t, shape1 = x$shape1, shape2 = x$shape2, log = log)
  }
}

#' Cumulative distribution function for a beta distribution.
#'
#' Returns a function that evaluates the beta CDF at given points.
#'
#' @param x A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(q, log.p = FALSE, ...)} returning the
#'   CDF (or log-CDF) at \code{q}.
#' @examples
#' x <- beta_dist(2, 5)
#' F <- cdf(x)
#' F(0.3)
#' F(0.5)
#' @importFrom stats pbeta
#' @export
cdf.beta_dist <- function(x, ...) {
  function(q, log.p = FALSE, ...) {
    pbeta(q, shape1 = x$shape1, shape2 = x$shape2, log.p = log.p)
  }
}

#' Inverse CDF (quantile function) for a beta distribution.
#'
#' Returns a function that computes quantiles of the beta distribution.
#'
#' @param x A \code{beta_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(p, lower.tail = TRUE, log.p = FALSE, ...)}
#'   returning the quantile at probability \code{p}.
#' @examples
#' x <- beta_dist(2, 5)
#' q <- inv_cdf(x)
#' q(0.5)
#' q(0.95)
#' @importFrom stats qbeta
#' @export
inv_cdf.beta_dist <- function(x, ...) {
  function(p, lower.tail = TRUE, log.p = FALSE, ...) {
    qbeta(p, shape1 = x$shape1, shape2 = x$shape2,
          lower.tail = lower.tail, log.p = log.p)
  }
}
