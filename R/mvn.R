#' Construct (multivariate or univariate) normal distribution object.
#'
#' @param mu mean vector
#' @param sigma variance-covariance matrix
#' @export
mvn <- function(
    mu = 0,
    sigma = diag(length(mu))) {

    stopifnot(is.numeric(mu),
              is.matrix(sigma),
              nrow(sigma) == ncol(sigma),
              nrow(sigma) == length(mu))

    structure(list(mu = mu, sigma = sigma),
              class = c("mvn"))
}

#' Retrieve the variance-covariance matrix (or scalar)
#' of a `mvn` object.

#' A variance-covariance matrix is a square matrix
#' giving the covariance between each pair of elements
#' of a given random vector. Importantly, its diagonal
#' contains variances of the elements
#'
#' @param object The `mvn` object to retrieve
#'               the variance-covariance matrix from
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `mvn` object
#' @importFrom stats vcov
#' @export
vcov.mvn <- function(object, ...) {
    object$sigma
}

#' Retrieve the mean of a `mvn` object.
#' @param x The `mvn` object to retrieve the mean from
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `mvn` object
#' @export
mean.mvn <- function(x, ...) {
    x$mu
}

#' Method for obtaining the parameters of a `mvn` object.
#'
#' @param x The object to obtain the parameters of
#' @param ... Additional arguments to pass (not used)
#' @return A named vector of parameters
#' @export
params.mvn <- function(x, ...) {
    c("mu" = x$mu, "sigma" = x$sigma)
}

#' Function to determine whether an object `x` is an `mvn` object.
#' @param x The object to test
#' @export
is_mvn <- function(x) {
    inherits(x, "mvn")
}

#' Method for sampling from a `mvn` (multivariate normal) object.
#'
#' @param x The `mvn` object to sample from
#' @param ... Additional arguments to pass to `rmvnorm` or `sample_mvn_region`
#' @importFrom mvtnorm rmvnorm
#' @export
sampler.mvn <- function(x, ...) {
    function(n = 1, mu = x$mu, sigma = x$sigma, p = 1) {
        if (p < 1) {
            sample_mvn_region(n, mu, sigma, p, ...)
        } else {
            rmvnorm(n, mu, sigma, ...)
        }
    }
}

#' Method for obtaining the pdf of an `mvn` object.
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function that computes the pdf of the mvn distribution.
#'         It accepts as input a parameter vector `x`, a mean vector `mu`,
#'         a variance-covariance matrix `sigma`, and a `log` argument
#'         determining whether to compute the log of the pdf. By default,
#'         `mu` and `sigma` are the mean and variance-covariance matrix
#'         of object `x`.
#' @importFrom mvtnorm dmvnorm
#' @export
pdf.mvn <- function(x, ...) {
    function(x, mu = x$mu, sigma = x$sigma, log = FALSE, ...) {
        dmvnorm(x = x, mean = mu, sigma = sigma, log = log, ...)
    }
}

#' Method for obtaining the support of a `mvn` object, where the support
#' is defined as values that have non-zero probability density.
#' 
#' @param x The `mvn` object to obtain the support of
#' @param ... Additional arguments to pass (not used)
#' @return A support-type object (see `support.R`), in this case an
#'         `interval` object for each component. 
#' @export
sup.mvn <- function(x, ...) {
    interval$new(lower = rep(-Inf, length(x$mu)),
                 upper = rep(Inf, length(x$mu)),
                 lower_closed = rep(FALSE, length(x$mu)),
                 upper_closed = rep(FALSE, length(x$mu)))
}

#' Generic method for obtaining the marginal distribution of an `mvn` object
#' `x` over components `indices`.
#' @param x The `mvn` object.
#' @param indices The indices of the marginal distribution to obtain.
#' @export
marginal.mvn <- function(x, indices) {
    if (length(indices) == 1) {
        return(normal(x$mu[indices], x$sigma[indices,indices]))
    }
    mu <- x$mu[indices]
    sigma <- x$sigma[indices, indices]
    mvn(mu, sigma)
}

#' Function for obtaining sample points for an `mvn` object that is within
#' the `p`-probability region.
#'
#' @param n the sample size
#' @param mu mean vector
#' @param sigma variance-covariance matrix
#' @param p the probability region
#' @param ... additional arguments to pass into `mahalanobis`
#' @importFrom stats qchisq mahalanobis
#' @export
sample_mvn_region <- function(n, mu, sigma, p = .95, ...) {

    stopifnot(p >= 0 && p <= 1, n > 0)

    k <- length(mu)
    crit <- qchisq(p, k)
    i <- 0L
    samp <- sampler(x)
    data <- matrix(nrow = n, ncol = k)
    while (i < n) {
        x <- samp(1)
        d <- mahalanobis(x, center = mu, cov = sigma, ...)
        if (d <= crit) {
            i <- i + 1L
            data[i, ] <- x
        }
    }
    data
}


#' Method for printing an `mvn` object.
#' @param x The object to print
#' @param ... Additional arguments to pass to `print`
#' @export
print.mvn <- function(x, ...) {
    cat("Multivariate normal distribution with mean:\n")
    print(x$mu)
    cat("and variance-covariance matrix:\n")
    print(x$sigma)
}


#' Method for obtaining the dimension of an `mvn` object.
#' @param x The object to obtain the dimension of
#' @return The dimension of the `mvn` object
#' @export
dim.mvn <- function(x) {
    length(x$mu)
}