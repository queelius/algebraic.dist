#' Construct a gamma distribution object.
#'
#' @param shape Shape parameter (positive scalar)
#' @param rate Rate parameter (positive scalar)
#' @return A `gamma_dist` object
#' @export
gamma_dist <- function(shape, rate) {
  if (!is.numeric(shape) || length(shape) != 1 || is.na(shape) || shape <= 0)
    stop("'shape' must be a positive scalar, got: ", deparse(shape))
  if (!is.numeric(rate) || length(rate) != 1 || is.na(rate) || rate <= 0)
    stop("'rate' must be a positive scalar, got: ", deparse(rate))
  structure(list(shape = shape, rate = rate),
            class = c("gamma_dist", "univariate_dist",
                      "continuous_dist", "dist"))
}

#' Test whether an object is a `gamma_dist`.
#' @param x The object to test
#' @return Logical; TRUE if `x` inherits from `gamma_dist`
#' @export
is_gamma_dist <- function(x) {
  inherits(x, "gamma_dist")
}

#' Method for obtaining the parameters of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @return A named numeric vector of parameters
#' @export
params.gamma_dist <- function(x) {
  c(shape = x$shape, rate = x$rate)
}

#' Retrieve the mean of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return The mean, `shape / rate`
#' @export
mean.gamma_dist <- function(x, ...) {
  x$shape / x$rate
}

#' Retrieve the variance of a `gamma_dist` object.
#'
#' @param object The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return The variance, `shape / rate^2`
#' @export
vcov.gamma_dist <- function(object, ...) {
  object$shape / object$rate^2
}

#' Retrieve the dimension of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @return 1 (univariate)
#' @export
dim.gamma_dist <- function(x) {
  1
}

#' Format a `gamma_dist` object as a character string.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A character string describing the distribution
#' @export
format.gamma_dist <- function(x, ...) {
  sprintf("Gamma distribution (shape = %g, rate = %g)", x$shape, x$rate)
}

#' Print method for `gamma_dist` objects.
#'
#' @param x The `gamma_dist` object to print
#' @param ... Additional arguments (not used)
#' @export
print.gamma_dist <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Method for sampling from a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object to sample from
#' @param ... Additional arguments (not used)
#' @return A function that generates `n` samples from the gamma distribution
#' @importFrom stats rgamma
#' @export
sampler.gamma_dist <- function(x, ...) {
  function(n = 1) {
    rgamma(n = n, shape = x$shape, rate = x$rate)
  }
}

#' Method for obtaining the density (pdf) of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A function that computes the pdf at point(s) `t`
#' @importFrom stats dgamma density
#' @export
density.gamma_dist <- function(x, ...) {
  function(t, log = FALSE) {
    dgamma(x = t, shape = x$shape, rate = x$rate, log = log)
  }
}

#' Method for obtaining the cdf of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A function that computes the cdf at point(s) `t`
#' @importFrom stats pgamma
#' @export
cdf.gamma_dist <- function(x, ...) {
  function(q, log.p = FALSE) {
    pgamma(q = q, shape = x$shape, rate = x$rate, log.p = log.p)
  }
}

#' Method for obtaining the inverse cdf (quantile function) of a `gamma_dist`
#' object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A function that computes the quantile at probability `p`
#' @importFrom stats qgamma
#' @export
inv_cdf.gamma_dist <- function(x, ...) {
  function(p, lower.tail = TRUE, log.p = FALSE) {
    qgamma(p = p, shape = x$shape, rate = x$rate,
           lower.tail = lower.tail, log.p = log.p)
  }
}

#' Support for gamma distribution, the positive real numbers (0, Inf).
#'
#' @param x The `gamma_dist` object
#' @return An `interval` object representing (0, Inf)
#' @export
sup.gamma_dist <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}

#' Method for obtaining the survival function of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A function that computes S(t) = P(X > t)
#' @export
surv.gamma_dist <- function(x, ...) {
  function(t, log.p = FALSE) {
    pgamma(q = t, shape = x$shape, rate = x$rate,
           lower.tail = FALSE, log.p = log.p)
  }
}

#' Method for obtaining the hazard function of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A function that computes h(t) = f(t) / S(t)
#' @export
hazard.gamma_dist <- function(x, ...) {
  function(t, log = FALSE) {
    log_f <- dgamma(t, shape = x$shape, rate = x$rate, log = TRUE)
    log_S <- pgamma(t, shape = x$shape, rate = x$rate,
                    lower.tail = FALSE, log.p = TRUE)
    log_h <- log_f - log_S
    if (log) log_h else exp(log_h)
  }
}
