#' Construct a Weibull distribution object.
#'
#' Creates an S3 object representing a Weibull distribution with the given
#' shape and scale parameters. The Weibull PDF is
#' \deqn{f(t) = (shape/scale)(t/scale)^{shape-1} \exp(-(t/scale)^{shape})}
#' for \eqn{t > 0}.
#'
#' @param shape Positive scalar shape parameter.
#' @param scale Positive scalar scale parameter.
#' @return A \code{weibull_dist} object with classes
#'   \code{c("weibull_dist", "univariate_dist", "continuous_dist", "dist")}.
#' @examples
#' x <- weibull_dist(shape = 2, scale = 3)
#' mean(x)
#' vcov(x)
#' format(x)
#' @export
weibull_dist <- function(shape, scale) {
  if (!is.numeric(shape) || length(shape) != 1 || is.na(shape) || shape <= 0)
    stop("'shape' must be a positive scalar, got: ", deparse(shape))
  if (!is.numeric(scale) || length(scale) != 1 || is.na(scale) || scale <= 0)
    stop("'scale' must be a positive scalar, got: ", deparse(scale))
  structure(list(shape = shape, scale = scale),
            class = c("weibull_dist", "univariate_dist",
                      "continuous_dist", "dist"))
}

#' Test whether an object is a \code{weibull_dist}.
#'
#' @param x The object to test.
#' @return \code{TRUE} if \code{x} inherits from \code{"weibull_dist"},
#'   \code{FALSE} otherwise.
#' @examples
#' is_weibull_dist(weibull_dist(2, 3))
#' is_weibull_dist(normal(0, 1))
#' @export
is_weibull_dist <- function(x) {
  inherits(x, "weibull_dist")
}

#' Retrieve the parameters of a \code{weibull_dist} object.
#'
#' @param x A \code{weibull_dist} object.
#' @return A named numeric vector with elements \code{shape} and \code{scale}.
#' @examples
#' params(weibull_dist(2, 3))
#' @export
params.weibull_dist <- function(x) {
  c("shape" = x$shape, "scale" = x$scale)
}

#' Mean of a Weibull distribution.
#'
#' Computes \eqn{scale \cdot \Gamma(1 + 1/shape)}.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return The mean of the distribution.
#' @examples
#' mean(weibull_dist(shape = 2, scale = 3))
#' @export
mean.weibull_dist <- function(x, ...) {
  x$scale * gamma(1 + 1 / x$shape)
}

#' Variance of a Weibull distribution.
#'
#' Computes \eqn{scale^2 (\Gamma(1 + 2/shape) - [\Gamma(1 + 1/shape)]^2)}.
#'
#' @param object A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return The variance (scalar).
#' @examples
#' vcov(weibull_dist(shape = 2, scale = 3))
#' @export
vcov.weibull_dist <- function(object, ...) {
  object$scale^2 * (gamma(1 + 2 / object$shape) -
                       gamma(1 + 1 / object$shape)^2)
}

#' Dimension of a Weibull distribution (always 1).
#'
#' @param x A \code{weibull_dist} object.
#' @return \code{1}.
#' @examples
#' dim(weibull_dist(2, 3))
#' @export
dim.weibull_dist <- function(x) {
  1
}

#' Format a \code{weibull_dist} object as a character string.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return A character string describing the distribution.
#' @examples
#' format(weibull_dist(2, 3))
#' @export
format.weibull_dist <- function(x, ...) {
  sprintf("Weibull distribution (shape = %g, scale = %g)", x$shape, x$scale)
}

#' Print a \code{weibull_dist} object.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return \code{x}, invisibly.
#' @examples
#' print(weibull_dist(2, 3))
#' @export
print.weibull_dist <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Sampler for a Weibull distribution.
#'
#' Returns a function that draws \code{n} independent samples from the
#' Weibull distribution.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(n = 1, ...)} returning a numeric vector
#'   of length \code{n}.
#' @examples
#' x <- weibull_dist(shape = 2, scale = 3)
#' s <- sampler(x)
#' set.seed(42)
#' s(5)
#' @importFrom stats rweibull
#' @export
sampler.weibull_dist <- function(x, ...) {
  function(n = 1, ...) {
    rweibull(n, shape = x$shape, scale = x$scale)
  }
}

#' Probability density function for a Weibull distribution.
#'
#' Returns a function that evaluates the Weibull PDF at given points.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log = FALSE, ...)} returning the
#'   density (or log-density) at \code{t}.
#' @examples
#' x <- weibull_dist(shape = 2, scale = 3)
#' f <- density(x)
#' f(1)
#' f(3)
#' @importFrom stats dweibull density
#' @export
density.weibull_dist <- function(x, ...) {
  function(t, log = FALSE, ...) {
    dweibull(t, shape = x$shape, scale = x$scale, log = log)
  }
}

#' Cumulative distribution function for a Weibull distribution.
#'
#' Returns a function that evaluates the Weibull CDF at given points.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(q, log.p = FALSE, ...)} returning the
#'   CDF (or log-CDF) at \code{q}.
#' @examples
#' x <- weibull_dist(shape = 2, scale = 3)
#' F <- cdf(x)
#' F(1)
#' F(3)
#' @importFrom stats pweibull
#' @export
cdf.weibull_dist <- function(x, ...) {
  function(q, log.p = FALSE, ...) {
    pweibull(q, shape = x$shape, scale = x$scale, log.p = log.p)
  }
}

#' Inverse CDF (quantile function) for a Weibull distribution.
#'
#' Returns a function that computes quantiles of the Weibull distribution.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(p, lower.tail = TRUE, log.p = FALSE, ...)}
#'   returning the quantile at probability \code{p}.
#' @examples
#' x <- weibull_dist(shape = 2, scale = 3)
#' q <- inv_cdf(x)
#' q(0.5)
#' q(0.95)
#' @importFrom stats qweibull
#' @export
inv_cdf.weibull_dist <- function(x, ...) {
  function(p, lower.tail = TRUE, log.p = FALSE, ...) {
    qweibull(p, shape = x$shape, scale = x$scale,
             lower.tail = lower.tail, log.p = log.p)
  }
}

#' Survival function for a Weibull distribution.
#'
#' Returns a function that computes \eqn{S(t) = P(X > t)} for the Weibull
#' distribution.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log.p = FALSE, ...)} returning the
#'   survival probability (or log-survival probability) at \code{t}.
#' @examples
#' x <- weibull_dist(shape = 2, scale = 3)
#' S <- surv(x)
#' S(1)
#' S(3)
#' @importFrom stats pweibull
#' @export
surv.weibull_dist <- function(x, ...) {
  function(t, log.p = FALSE, ...) {
    pweibull(t, shape = x$shape, scale = x$scale,
             lower.tail = FALSE, log.p = log.p)
  }
}

#' Hazard function for a Weibull distribution.
#'
#' Returns a function that evaluates the Weibull hazard rate
#' \eqn{h(t) = (shape/scale)(t/scale)^{shape-1}} for \eqn{t > 0}.
#'
#' @param x A \code{weibull_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log = FALSE)} returning the hazard
#'   (or log-hazard) at \code{t}.
#' @examples
#' x <- weibull_dist(shape = 2, scale = 3)
#' h <- hazard(x)
#' h(1)
#' h(3)
#' @export
hazard.weibull_dist <- function(x, ...) {
  function(t, log = FALSE) {
    h <- ifelse(t <= 0, 0, (x$shape / x$scale) * (t / x$scale)^(x$shape - 1))
    if (log) {
      ifelse(t <= 0, -Inf, log(h))
    } else {
      h
    }
  }
}

#' Support of a Weibull distribution.
#'
#' The Weibull distribution is supported on \eqn{(0, \infty)}.
#'
#' @param x A \code{weibull_dist} object.
#' @return An \code{interval} object representing \eqn{(0, \infty)}.
#' @examples
#' sup(weibull_dist(2, 3))
#' @export
sup.weibull_dist <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}
