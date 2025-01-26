#' Construct univariate normal distribution object.
#'
#' @param mu mean
#' @param var variance
#' @export
normal <- function(mu = 0, var = 1) {
    stopifnot(is.numeric(mu), is.numeric(var))
    structure(list(mu = mu, var = var),
              class = c("normal", "univariate_dist",
                        "continuous_dist", "dist"))
}

#' Method for printing a `normal` object.
#' @param x The object to print
#' @param ... Additional arguments to pass (not used)
#' @export
print.normal <- function(x, ...) {
    cat("Normal distribution with mean", x$mu, "and variance", x$var, "\n")
}

#' Retrieve the variance-covariance matrix (or scalar)
#' of a `normal` object.
#'
#' @param object The `normal` object to retrieve the variance-covariance matrix from
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `normal` object
#' @export
vcov.normal <- function(object, ...) {
    object$var
}

#' Retrieve the mean of a `normal` object.
#' @param x The `normal` object to retrieve the mean from
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `normal` object
#' @export
mean.normal <- function(x, ...) {
    x$mu
}

#' Method for obtaining the parameters of a `normal` object.
#'
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @export
params.normal <- function(x) {
    c("mu" = x$mu, "var" = x$var)
}

#' Function to determine whether an object `x` is an `normal` object.
#' @param x The object to test
#' @export
is_normal <- function(x) {
    inherits(x, "normal")
}

#' Method for sampling from a `normal` object.
#'
#' @param x The `normal` object to sample from
#' @param ... Additional arguments to pass (not used)
#' @return A function that samples from the normal distribution. As input,
#'         it accepts a sample size `n`, a numeric `mu`, and a variance
#'         numeric `var`. By default, `mu` and `var` are the mean and
#'         variance of object `x`.
#' @importFrom stats rnorm
#' @export
sampler.normal <- function(x, ...) {
    function(n = 1, mu = x$mu, var = x$var, ...) {
        rnorm(n, mu, var, ...)
    }
}

#' Method for obtaining the pdf of an `normal` object.
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function that computes the pdf of the normal distribution.
#'         It accepts as input a parameter vector `x`, a mean vector `mu`,
#'         a variance-covariance matrix `var`, and a `log` argument
#'         determining whether to compute the log of the pdf. By default,
#'         `mu` and `var` are the mean and variance of object `x`.
#' @importFrom stats dnorm density
#' @export
density.normal <- function(x, ...) {
    function(t, mu = x$mu, var = x$var, log = FALSE, ...) {
        dnorm(x = t, mean = mu, sd = sqrt(var), log = log, ...)
    }
}

#' Method for obtaining the cdf of an `normal` object.
#' @param x The object to obtain the cdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function that computes the cdf of the normal distribution.
#'         It accepts as input a parameter vector `q`, a mean vector `mu`,
#'         a variance `var`, and a `log` argument
#'         determining whether to compute the log of the cdf. By default,
#'         `mu` and `var` are the mean and variance of object `x` and `log`
#'         is `FALSE`. Finally, it accepts additional arguments `...` to
#'         pass to the `pnorm` function.
#' @importFrom stats pnorm
#' @export
cdf.normal <- function(x, ...) {
    function(q, mu = x$mu, var = x$var, log.p = FALSE, ...) {
        pnorm(q = q, mean = mu, sd = sqrt(var), log.p = log.p, ...)
    }
}

#' Method for obtaining the support of a `normal` object, where the support
#' is defined as values that have non-zero probability density.
#' 
#' @param x The `normal` object to obtain the support of
#' @return A support-type object (see `support.R`), in this case an
#'         `interval` object for each component. 
#' @export
sup.normal <- function(x) {
    interval$new()
}

#' Method for obtaining the dimension of a `normal` object.
#' @param x The `normal` object to obtain the dimension of
#' @return The dimension of the `normal` object
#' @export
dim.normal <- function(x) {
    1
}

#' Method for obtaining the inverse cdf of an `normal` object.
#' @param x The object to obtain the inverse cdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function that computes the inverse cdf of the normal distribution.
#' @importFrom stats qnorm
#' @export
inv_cdf.normal <- function(x, ...) {
    function(p, mu = x$mu, var = x$var, log.p = FALSE, ...) {
        qnorm(p, mu, sqrt(var), ...)
    }
}

#' Operator to take the negative of a `normal` object.
#' 
#' @param e The `normal` object to take the negative of
#' @return A new `normal` object with mean equal to the negative of the mean
#'         of the input object and variance equal to the variance of the input
#'         object.
#' @export
`-.normal` <- function(e) {
    normal(mu = -mean(e), var = vcov(e))
}

