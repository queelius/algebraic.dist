#' Construct a multivariate or univariate normal distribution object.
#'
#' This function constructs an object representing a normal distribution.
#' If the length of the mean vector `mu` is 1, it creates a univariate 
#' normal distribution. Otherwise, it creates a multivariate normal distribution.
#'
#' @param mu A numeric vector specifying the means of the distribution. 
#'           If `mu` has length 1, a univariate normal distribution is created. 
#'           If `mu` has length > 1, a multivariate normal distribution is created.
#' @param sigma A numeric matrix specifying the variance-covariance matrix of the 
#'              distribution. It must be a square matrix with the same number of 
#'              rows and columns as the length of `mu`. Default is the identity 
#'              matrix of size equal to the length of `mu`.
#' @return If `mu` has length 1, it returns a `normal` object. If `mu` has length 
#'         > 1, it returns an `mvn` object. Both types of objects contain `mu` 
#'         and `sigma` as their properties.
#' @export
mvn <- function(
    mu,
    sigma = diag(length(mu))) {

    stopifnot(is.numeric(mu),
              is.matrix(sigma),
              nrow(sigma) == ncol(sigma),
              nrow(sigma) == length(mu))

    if (length(mu) == 1) {
        return(normal(mu, sigma))
    } else {
        structure(list(mu = mu, sigma = sigma),
            class = c("mvn", "multivariate_dist",
                      "continuous_dist", "dist"))
    }
}

#' Retrieve the variance-covariance matrix of an `mvn` object.
#' @param object The `mvn` object to retrieve the variance-covariance matrix of
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `mvn` object
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
#' @return A named vector of parameters
#' @export
params.mvn <- function(x) {
    c("mu" = x$mu, "sigma" = x$sigma)
}

#' Function to determine whether an object `x` is an `mvn` object.
#' @param x The object to test
#' @export
is_mvn <- function(x) {
    inherits(x, "mvn")
}

#' Function generator for sampling from a `mvn` (multivariate normal) object.
#'
#' @param x The `mvn` object to sample from
#' @param ... Additional arguments to pass to the generated function that will
#' be fixed during all calls.
#' @return A function that samples from the `mvn` distribution. It accepts as
#'         input:
#'          - `n`: number of samples to generate. Defaults to 1.
#'          - `mu`: a vector denoting the population mean. Defaults to the mean 
#'            of `x` (an `mvn` object)
#'          - `sigma`: a matrix denoting the covariance of observations. 
#'            Defaults to the variance-covariance of `x`.
#'          - `p`: probability region to sample from. Defaults to 1, which
#'            corresponds to the entire distribution.
#'            `sample_mvn_region` method. It's used when `p` is less than 1.
#'          - `...`: any additional parameters to pass to `rmvnorm` or 
#'            `sample_mvn_region` which can be different during each call.
#' @importFrom mvtnorm rmvnorm
#' @export
sampler.mvn <- function(x, ...) {
    fixed_args <- list(...)
    function(n = 1, mu = x$mu, sigma = x$sigma, p = 1, ...) {
        if (p < 1) {
            do.call(sample_mvn_region, c(list(n = n, mean = mu, sigma = sigma,
                p = p), fixed_args, list(...)))
        } else {
            do.call(rmvnorm, c(list(n = n, mean = mu, sigma = sigma),
                fixed_args, list(...)))
        }
    }
}

#' Function generator for obtaining the pdf of an `mvn` object (multivariate
#' normal).
#' 
#' @param x The `mvn` (S3) object to obtain the pdf (density) of
#' @param ... Additional arguments to pass into the generated function.
#' @return A function that computes the pdf of the `mvn` distribution.
#'         It accepts as input:
#'          - `obs`: vector or matrix of quantiles. when x is a matrix, each row
#'            is taken to be a quantile and columns correspond to the number of
#'            dimensions, p.
#'          - `mu`: a a vector denoting the population mean. Defaults to the
#'            mean of `x` (an `mvn` object)
#'          - `sigma`: a matrix denoting the variance-covariance of
#'            observations. Defaults to the variance-covariance of `x`.
#'          - `log`: logical, determines whether to compute the log of the pdf.
#'            Defaults to `FALSE`.
#'          - `...`: any additional parameters to pass to `dmvnorm`.
#' @importFrom mvtnorm dmvnorm
#' @export
density.mvn <- function(x, ...) {
    fixed_args <- list(...)
    function(obs, mu = x$mu, sigma = x$sigma, log = FALSE, ...) {
        do.call(dmvnorm, c(list(x = obs, mean = mu, sigma = sigma, log = log),
            fixed_args, list(...)))
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
    if (length(indices) == 0) {
        stop("indices must be non-empty")
    }
    else if (length(indices) == 1) {
        return(normal(x$mu[indices], x$sigma[indices,indices]))
    }
    else if (any(indices < 0) || any(indices > dim(x))) {
        stop("indices must be in [1, dim(x)]")
    }
    mu <- x$mu[indices]
    sigma <- x$sigma[indices, indices]
    mvn(mu, sigma)
}

#' Function for obtaining sample points for an `mvn` object that is within
#' the `p`-probability region. That is, it samples from the smallest region of
#' the distribution that contains `p` probability mass. This is done by first
#' sampling from the entire distribution, then rejecting samples that are not
#' in the probability region (using the statistical distance `mahalanobis`
#' from `mu`).
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

#' Method for obtaining the CDF of a `mvn` object.
#' @param x The object to obtain the CDF of
#' @param ... Additional arguments to pass (not used)
#' @importFrom mvtnorm pmvnorm
#' @export
cdf.mvn <- function(x, ...) {
    function(q, mu = x$mu, sigma = x$sigma, lower.tail = TRUE,
             log.p = FALSE, ...) {
        p <- pmvnorm(mean = mu, sigma = sigma, upper = q, ...)
        # check to see if `p` numeric has an attr `error` set to non-zero
        if (attr(p, "error") != 0) {
            warning(paste0("error (", attr(p, "error"),"): ", attr(p, "msg")))
        }
        attr(p, "error") <- NULL
        attr(p, "msg") <- NULL
        if (log.p) {
            return(log(p))
        } else {
            return(p)
        }
    }
}

#' Computes the distribution of `g(x)` where `x` is an `mvn` object.
#'
#' By the invariance property, if `x` is an `mvn` object,
#' then under the right conditions, asymptotically, `g(x)` is an MVN
#' distributed,
#'     g(x) ~ normal(g(mean(x)), sigma)
#' where `sigma` is the variance-covariance of `g(x)`
#'
#' @param x The `mvn` object to apply `g` to
#' @param g The function to apply to `x`
#' @param n number of samples to take to estimate distribution of `g(x)` if
#'         `method` is `mc` or `empirical`. Defaults to 10000.
#' @param ... additional arguments to pass into the `g` function.
#' @importFrom stats vcov
#' @export
rmap.mvn <- function(x, g, n = 10000L, ...) {
    D <- rmap(empirical_dist(sampler(x)(n)), g, ...)
    mvn(mu = mean(D), sigma = vcov(D))
}
