#' Construct univariate normal distribution object.
#'
#' @param mu mean
#' @param var variance
#' @export
normal <- function(
    mu = 0,
    var = 1) {

    stopifnot(is.numeric(mu),
              is.numeric(var))

    structure(list(mu = mu, var = var),
              class = c("normal", "univariate_dist", "dist"))
}

#' Retrieve the variance-covariance matrix (or scalar)
#' of a `normal` object.

#' A variance-covariance matrix is a square matrix
#' giving the covariance between each pair of elements
#' of a given random vector. Importantly, its diagonal
#' contains variances of the elements
#'
#' @param object The `normal` object to retrieve
#'               the variance-covariance matrix from
#' @return The variance-covariance matrix of the `normal` object
#' @importFrom stats vcov
#' @export
vcov.normal <- function(object) {
    object$var
}

#' Retrieve the mean of a `normal` object.
#' @param object The `normal` object to retrieve the mean from
#' @return The mean of the `normal` object
#' @export
mean.normal <- function(object) {
    object$mu
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
#' @return A function that samples from the normal distribution. As input,
#'         it accepts a sample size `n`, a numeric `mu`, and a variance
#'         numeric `var`. By default, `mu` and `var` are the mean and
#'         variance of object `x`.
#' @importFrom stats rnorm
#' @export
sampler.normal <- function(x) {
    function(n = 1, mu = x$mu, var = x$var) {
        rnorm(n, mu, var)
    }
}

#' Method for obtaining the pdf of an `normal` object.
#' @param x The object to obtain the pdf of
#' @return A function that computes the pdf of the normal distribution.
#'         It accepts as input a parameter vector `x`, a mean vector `mu`,
#'         a variance-covariance matrix `var`, and a `log` argument
#'         determining whether to compute the log of the pdf. By default,
#'         `mu` and `var` are the mean and variance of object `x`.
#' @importFrom stats dnorm
#' @export
pdf.normal <- function(x) {
    function(x, mu = x$mu, var = x$var, log = FALSE) {
        dnorm(x = x, mean = mu, sd = sqrt(var), log = log)
    }
}

#' Method for obtaining the cdf of an `normal` object.
#' @param x The object to obtain the cdf of
#' @return A function that computes the cdf of the normal distribution.
#'         It accepts as input a parameter vector `x`, a mean vector `mu`,
#'         a variance-covariance matrix `var`, and a `log` argument
#'         determining whether to compute the log of the cdf. By default,
#'         `mu` and `var` are the mean and variance of object `x`.
#' @importFrom stats pnorm
#' @export
cdf.normal <- function(x) {
    function(q, mu = x$mu, var = x$var, log = FALSE) {
        pnorm(q = q, mean = mu, sd = sqrt(var), log = log)
    }
}

#' Method for obtaining the survival function of an `normal` object.
#' @param x The object to obtain the survival function of
#' @return A function that computes the survival of the normal distribution.
#'         It accepts as input a parameter vector `x`, a mean vector `mu`,
#'         a variance-covariance matrix `sigma`, and a `log` argument
#'         determining whether to compute the log of the cdf. By default,
#'         `mu` and `sigma` are the mean and variance of object `x`.
#' @importFrom stats pnorm
#' @export
surv.normal <- function(x) {
    function(x, mu = x$mu, var = x$var, log = FALSE) {
        pnorm(x = x, mean = mu, sd = sqrt(var), log = log, lower.tail = FALSE)
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
