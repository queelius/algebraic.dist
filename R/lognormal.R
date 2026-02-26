#' Construct a log-normal distribution object.
#'
#' Creates an S3 object representing a log-normal distribution with the given
#' \code{meanlog} and \code{sdlog} parameters. The log-normal PDF is
#' \deqn{f(t) = \frac{1}{t \cdot sdlog \sqrt{2\pi}}
#'   \exp\!\left(-\frac{(\log t - meanlog)^2}{2 \cdot sdlog^2}\right)}
#' for \eqn{t > 0}.
#'
#' @param meanlog Mean of the distribution on the log scale (default 0).
#' @param sdlog Standard deviation on the log scale (default 1), must be
#'   positive.
#' @return A \code{lognormal} object with classes
#'   \code{c("lognormal", "univariate_dist", "continuous_dist", "dist")}.
#' @export
lognormal <- function(meanlog = 0, sdlog = 1) {
  if (!is.numeric(meanlog) || length(meanlog) != 1 || is.na(meanlog))
    stop("'meanlog' must be a numeric scalar, got: ", deparse(meanlog))
  if (!is.numeric(sdlog) || length(sdlog) != 1 || is.na(sdlog) || sdlog <= 0)
    stop("'sdlog' must be a positive scalar, got: ", deparse(sdlog))
  structure(list(meanlog = meanlog, sdlog = sdlog),
            class = c("lognormal", "univariate_dist",
                      "continuous_dist", "dist"))
}

#' Test whether an object is a \code{lognormal}.
#'
#' @param x The object to test.
#' @return \code{TRUE} if \code{x} inherits from \code{"lognormal"},
#'   \code{FALSE} otherwise.
#' @export
is_lognormal <- function(x) {
  inherits(x, "lognormal")
}

#' Retrieve the parameters of a \code{lognormal} object.
#'
#' @param x A \code{lognormal} object.
#' @return A named numeric vector with elements \code{meanlog} and \code{sdlog}.
#' @export
params.lognormal <- function(x) {
  c(meanlog = x$meanlog, sdlog = x$sdlog)
}

#' Mean of a log-normal distribution.
#'
#' Computes \eqn{\exp(meanlog + sdlog^2 / 2)}.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return The mean of the distribution.
#' @export
mean.lognormal <- function(x, ...) {
  exp(x$meanlog + x$sdlog^2 / 2)
}

#' Variance of a log-normal distribution.
#'
#' Computes \eqn{(\exp(sdlog^2) - 1) \exp(2 \cdot meanlog + sdlog^2)}.
#'
#' @param object A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return The variance (scalar).
#' @export
vcov.lognormal <- function(object, ...) {
  (exp(object$sdlog^2) - 1) * exp(2 * object$meanlog + object$sdlog^2)
}

#' Dimension of a log-normal distribution (always 1).
#'
#' @param x A \code{lognormal} object.
#' @return \code{1}.
#' @export
dim.lognormal <- function(x) {
  1
}

#' Format a \code{lognormal} object as a character string.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return A character string describing the distribution.
#' @export
format.lognormal <- function(x, ...) {
  sprintf("Log-normal distribution (meanlog = %g, sdlog = %g)",
          x$meanlog, x$sdlog)
}

#' Print a \code{lognormal} object.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return \code{x}, invisibly.
#' @export
print.lognormal <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Support of a log-normal distribution.
#'
#' The log-normal distribution is supported on \eqn{(0, \infty)}.
#'
#' @param x A \code{lognormal} object.
#' @return An \code{interval} object representing \eqn{(0, \infty)}.
#' @export
sup.lognormal <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}

#' Sampler for a log-normal distribution.
#'
#' Returns a function that draws \code{n} independent samples from the
#' log-normal distribution.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(n = 1, ...)} returning a numeric vector
#'   of length \code{n}.
#' @importFrom stats rlnorm
#' @export
sampler.lognormal <- function(x, ...) {
  function(n = 1, ...) {
    rlnorm(n, meanlog = x$meanlog, sdlog = x$sdlog)
  }
}

#' Probability density function for a log-normal distribution.
#'
#' Returns a function that evaluates the log-normal PDF at given points.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log = FALSE, ...)} returning the
#'   density (or log-density) at \code{t}.
#' @importFrom stats dlnorm density
#' @export
density.lognormal <- function(x, ...) {
  function(t, log = FALSE, ...) {
    dlnorm(t, meanlog = x$meanlog, sdlog = x$sdlog, log = log)
  }
}

#' Cumulative distribution function for a log-normal distribution.
#'
#' Returns a function that evaluates the log-normal CDF at given points.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(q, log.p = FALSE, ...)} returning the
#'   CDF (or log-CDF) at \code{q}.
#' @importFrom stats plnorm
#' @export
cdf.lognormal <- function(x, ...) {
  function(q, log.p = FALSE, ...) {
    plnorm(q, meanlog = x$meanlog, sdlog = x$sdlog, log.p = log.p)
  }
}

#' Inverse CDF (quantile function) for a log-normal distribution.
#'
#' Returns a function that computes quantiles of the log-normal distribution.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(p, lower.tail = TRUE, log.p = FALSE, ...)}
#'   returning the quantile at probability \code{p}.
#' @importFrom stats qlnorm
#' @export
inv_cdf.lognormal <- function(x, ...) {
  function(p, lower.tail = TRUE, log.p = FALSE, ...) {
    qlnorm(p, meanlog = x$meanlog, sdlog = x$sdlog,
           lower.tail = lower.tail, log.p = log.p)
  }
}

#' Survival function for a log-normal distribution.
#'
#' Returns a function that computes \eqn{S(t) = P(X > t)} for the log-normal
#' distribution.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log.p = FALSE, ...)} returning the
#'   survival probability (or log-survival probability) at \code{t}.
#' @importFrom stats plnorm
#' @export
surv.lognormal <- function(x, ...) {
  function(t, log.p = FALSE, ...) {
    plnorm(t, meanlog = x$meanlog, sdlog = x$sdlog,
           lower.tail = FALSE, log.p = log.p)
  }
}

#' Hazard function for a log-normal distribution.
#'
#' Returns a function that evaluates the log-normal hazard rate
#' \eqn{h(t) = f(t) / S(t)} for \eqn{t > 0}.
#'
#' @param x A \code{lognormal} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log = FALSE)} returning the hazard
#'   (or log-hazard) at \code{t}.
#' @export
hazard.lognormal <- function(x, ...) {
  f <- density(x)
  S <- surv(x)
  function(t, log = FALSE) {
    h <- f(t) / S(t)
    if (log) log(h) else h
  }
}
