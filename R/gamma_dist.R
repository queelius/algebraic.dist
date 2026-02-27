#' Construct a gamma distribution object.
#'
#' @param shape Shape parameter (positive scalar)
#' @param rate Rate parameter (positive scalar)
#' @return A `gamma_dist` object
#' @examples
#' x <- gamma_dist(shape = 2, rate = 1)
#' mean(x)
#' vcov(x)
#' format(x)
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
#' @examples
#' is_gamma_dist(gamma_dist(2, 1))
#' is_gamma_dist(normal(0, 1))
#' @export
is_gamma_dist <- function(x) {
  inherits(x, "gamma_dist")
}

#' Method for obtaining the parameters of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @return A named numeric vector of parameters
#' @examples
#' params(gamma_dist(2, 1))
#' @export
params.gamma_dist <- function(x) {
  c(shape = x$shape, rate = x$rate)
}

#' Retrieve the mean of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return The mean, `shape / rate`
#' @examples
#' mean(gamma_dist(shape = 3, rate = 2))
#' @export
mean.gamma_dist <- function(x, ...) {
  x$shape / x$rate
}

#' Retrieve the variance of a `gamma_dist` object.
#'
#' @param object The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return The variance, `shape / rate^2`
#' @examples
#' vcov(gamma_dist(shape = 3, rate = 2))
#' @export
vcov.gamma_dist <- function(object, ...) {
  object$shape / object$rate^2
}

#' Retrieve the dimension of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @return 1 (univariate)
#' @examples
#' dim(gamma_dist(2, 1))
#' @export
dim.gamma_dist <- function(x) {
  1
}

#' Format a `gamma_dist` object as a character string.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A character string describing the distribution
#' @examples
#' format(gamma_dist(2, 1))
#' @export
format.gamma_dist <- function(x, ...) {
  sprintf("Gamma distribution (shape = %g, rate = %g)", x$shape, x$rate)
}

#' Print method for `gamma_dist` objects.
#'
#' @param x The `gamma_dist` object to print
#' @param ... Additional arguments (not used)
#' @return \code{x}, invisibly.
#' @examples
#' print(gamma_dist(2, 1))
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
#' @examples
#' x <- gamma_dist(shape = 2, rate = 1)
#' s <- sampler(x)
#' set.seed(42)
#' s(5)
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
#' @examples
#' x <- gamma_dist(shape = 2, rate = 1)
#' f <- density(x)
#' f(1)
#' f(2)
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
#' @examples
#' x <- gamma_dist(shape = 2, rate = 1)
#' F <- cdf(x)
#' F(1)
#' F(2)
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
#' @examples
#' x <- gamma_dist(shape = 2, rate = 1)
#' q <- inv_cdf(x)
#' q(0.5)
#' q(0.95)
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
#' @examples
#' sup(gamma_dist(2, 1))
#' @export
sup.gamma_dist <- function(x) {
  interval$new(lower = 0, lower_closed = FALSE)
}

#' Method for obtaining the survival function of a `gamma_dist` object.
#'
#' @param x The `gamma_dist` object
#' @param ... Additional arguments (not used)
#' @return A function that computes S(t) = P(X > t)
#' @examples
#' x <- gamma_dist(shape = 2, rate = 1)
#' S <- surv(x)
#' S(1)
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
#' @examples
#' x <- gamma_dist(shape = 2, rate = 1)
#' h <- hazard(x)
#' h(1)
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
