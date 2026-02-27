#' Construct a uniform distribution object.
#'
#' Creates an S3 object representing a continuous uniform distribution on the
#' interval \eqn{[min, max]}.  The PDF is \eqn{f(x) = 1/(max - min)} for
#' \eqn{min \le x \le max}.
#'
#' @param min Lower bound of the distribution (default 0).
#' @param max Upper bound of the distribution (default 1).
#' @return A \code{uniform_dist} object with classes
#'   \code{c("uniform_dist", "univariate_dist", "continuous_dist", "dist")}.
#' @examples
#' x <- uniform_dist(min = 0, max = 10)
#' mean(x)
#' vcov(x)
#' format(x)
#' @export
uniform_dist <- function(min = 0, max = 1) {
  if (!is.numeric(min) || length(min) != 1 || is.na(min))
    stop("'min' must be a numeric scalar, got: ", deparse(min))
  if (!is.numeric(max) || length(max) != 1 || is.na(max))
    stop("'max' must be a numeric scalar, got: ", deparse(max))
  if (min >= max)
    stop("'min' must be less than 'max', got min = ", min, ", max = ", max)
  structure(list(min = min, max = max),
            class = c("uniform_dist", "univariate_dist",
                      "continuous_dist", "dist"))
}

#' Test whether an object is a \code{uniform_dist}.
#'
#' @param x The object to test.
#' @return \code{TRUE} if \code{x} inherits from \code{"uniform_dist"},
#'   \code{FALSE} otherwise.
#' @examples
#' is_uniform_dist(uniform_dist(0, 1))
#' is_uniform_dist(normal(0, 1))
#' @export
is_uniform_dist <- function(x) {
  inherits(x, "uniform_dist")
}

#' Retrieve the parameters of a \code{uniform_dist} object.
#'
#' @param x A \code{uniform_dist} object.
#' @return A named numeric vector with elements \code{min} and \code{max}.
#' @examples
#' params(uniform_dist(0, 10))
#' @export
params.uniform_dist <- function(x) {
  c("min" = x$min, "max" = x$max)
}

#' Mean of a uniform distribution.
#'
#' Computes \eqn{(min + max) / 2}.
#'
#' @param x A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return The mean of the distribution.
#' @examples
#' mean(uniform_dist(0, 10))
#' @export
mean.uniform_dist <- function(x, ...) {
  (x$min + x$max) / 2
}

#' Variance of a uniform distribution.
#'
#' Computes \eqn{(max - min)^2 / 12}.
#'
#' @param object A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return The variance (scalar).
#' @examples
#' vcov(uniform_dist(0, 10))
#' @export
vcov.uniform_dist <- function(object, ...) {
  (object$max - object$min)^2 / 12
}

#' Dimension of a uniform distribution (always 1).
#'
#' @param x A \code{uniform_dist} object.
#' @return \code{1}.
#' @examples
#' dim(uniform_dist(0, 1))
#' @export
dim.uniform_dist <- function(x) {
  1
}

#' Format a \code{uniform_dist} object as a character string.
#'
#' @param x A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return A character string describing the distribution.
#' @examples
#' format(uniform_dist(0, 10))
#' @export
format.uniform_dist <- function(x, ...) {
  sprintf("Uniform distribution (min = %g, max = %g)", x$min, x$max)
}

#' Print a \code{uniform_dist} object.
#'
#' @param x A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return \code{x}, invisibly.
#' @examples
#' print(uniform_dist(0, 10))
#' @export
print.uniform_dist <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Support of a uniform distribution.
#'
#' The uniform distribution is supported on \eqn{[min, max]}.
#'
#' @param x A \code{uniform_dist} object.
#' @return An \code{interval} object representing \eqn{[min, max]}.
#' @examples
#' sup(uniform_dist(0, 10))
#' @export
sup.uniform_dist <- function(x) {
  interval$new(lower = x$min, upper = x$max,
               lower_closed = TRUE, upper_closed = TRUE)
}

#' Sampler for a uniform distribution.
#'
#' Returns a function that draws \code{n} independent samples from the
#' uniform distribution.
#'
#' @param x A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(n = 1, ...)} returning a numeric vector
#'   of length \code{n}.
#' @examples
#' x <- uniform_dist(0, 10)
#' s <- sampler(x)
#' set.seed(42)
#' s(5)
#' @importFrom stats runif
#' @export
sampler.uniform_dist <- function(x, ...) {
  function(n = 1, ...) {
    runif(n, min = x$min, max = x$max)
  }
}

#' Probability density function for a uniform distribution.
#'
#' Returns a function that evaluates the uniform PDF at given points.
#'
#' @param x A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log = FALSE, ...)} returning the
#'   density (or log-density) at \code{t}.
#' @examples
#' x <- uniform_dist(0, 10)
#' f <- density(x)
#' f(5)
#' f(15)
#' @importFrom stats dunif density
#' @export
density.uniform_dist <- function(x, ...) {
  function(t, log = FALSE, ...) {
    dunif(t, min = x$min, max = x$max, log = log)
  }
}

#' Cumulative distribution function for a uniform distribution.
#'
#' Returns a function that evaluates the uniform CDF at given points.
#'
#' @param x A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(q, log.p = FALSE, ...)} returning the
#'   CDF (or log-CDF) at \code{q}.
#' @examples
#' x <- uniform_dist(0, 10)
#' F <- cdf(x)
#' F(5)
#' F(10)
#' @importFrom stats punif
#' @export
cdf.uniform_dist <- function(x, ...) {
  function(q, log.p = FALSE, ...) {
    punif(q, min = x$min, max = x$max, log.p = log.p)
  }
}

#' Inverse CDF (quantile function) for a uniform distribution.
#'
#' Returns a function that computes quantiles of the uniform distribution.
#'
#' @param x A \code{uniform_dist} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(p, lower.tail = TRUE, log.p = FALSE, ...)}
#'   returning the quantile at probability \code{p}.
#' @examples
#' x <- uniform_dist(0, 10)
#' q <- inv_cdf(x)
#' q(0.5)
#' q(0.9)
#' @importFrom stats qunif
#' @export
inv_cdf.uniform_dist <- function(x, ...) {
  function(p, lower.tail = TRUE, log.p = FALSE, ...) {
    qunif(p, min = x$min, max = x$max,
           lower.tail = lower.tail, log.p = log.p)
  }
}
