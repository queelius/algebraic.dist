#' Construct univariate normal distribution object.
#'
#' @param mu mean
#' @param var variance
#' @return A \code{normal} distribution object.
#' @examples
#' x <- normal(mu = 0, var = 1)
#' mean(x)
#' vcov(x)
#' format(x)
#' @export
normal <- function(mu = 0, var = 1) {
    if (!is.numeric(mu) || length(mu) != 1 || is.na(mu))
        stop("'mu' must be a numeric scalar, got: ", deparse(mu))
    if (!is.numeric(var) || length(var) != 1 || is.na(var) || var < 0)
        stop("'var' must be a non-negative numeric scalar, got: ", deparse(var))
    structure(list(mu = mu, var = var),
              class = c("normal", "univariate_dist",
                        "continuous_dist", "dist"))
}

#' Format method for `normal` objects.
#' @param x The object to format
#' @param ... Additional arguments (not used)
#' @return A character string
#' @examples
#' x <- normal(2, 3)
#' format(x)
#' @export
format.normal <- function(x, ...) {
    sprintf("Normal distribution (mu = %g, var = %g)", x$mu, x$var)
}

#' Print method for `normal` objects.
#' @param x The object to print
#' @param ... Additional arguments to pass (not used)
#' @return \code{x}, invisibly.
#' @examples
#' x <- normal(2, 3)
#' print(x)
#' @export
print.normal <- function(x, ...) {
    cat(format(x), "\n")
    invisible(x)
}

#' Retrieve the variance-covariance matrix (or scalar)
#' of a `normal` object.
#'
#' @param object The `normal` object to retrieve the variance-covariance matrix from
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `normal` object
#' @examples
#' x <- normal(0, 4)
#' vcov(x)
#' @export
vcov.normal <- function(object, ...) {
    object$var
}

#' Retrieve the mean of a `normal` object.
#' @param x The `normal` object to retrieve the mean from
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `normal` object
#' @examples
#' x <- normal(5, 2)
#' mean(x)
#' @export
mean.normal <- function(x, ...) {
    x$mu
}

#' Method for obtaining the parameters of a `normal` object.
#'
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @examples
#' x <- normal(3, 2)
#' params(x)
#' @export
params.normal <- function(x) {
    c("mu" = x$mu, "var" = x$var)
}

#' Function to determine whether an object `x` is an `normal` object.
#' @param x The object to test
#' @return Logical; \code{TRUE} if \code{x} is a \code{normal}.
#' @examples
#' is_normal(normal(0, 1))
#' is_normal(exponential(1))
#' @export
is_normal <- function(x) {
    inherits(x, "normal")
}

#' Method for sampling from a `normal` object.
#'
#' @param x The `normal` object to sample from
#' @param ... Additional arguments to pass (not used)
#' @return A function \code{function(n = 1, ...)} that draws \code{n}
#'   samples from the normal distribution.
#' @examples
#' x <- normal(0, 1)
#' s <- sampler(x)
#' set.seed(42)
#' s(5)
#' @importFrom stats rnorm
#' @export
sampler.normal <- function(x, ...) {
    mu <- x$mu
    sd <- sqrt(x$var)
    function(n = 1, ...) {
        rnorm(n, mean = mu, sd = sd)
    }
}

#' Method for obtaining the pdf of an `normal` object.
#' @param x The object to obtain the pdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function \code{function(t, log = FALSE, ...)} that computes the
#'   pdf (or log-pdf) of the normal distribution at \code{t}.
#' @examples
#' x <- normal(0, 1)
#' f <- density(x)
#' f(0)
#' f(1)
#' @importFrom stats dnorm density
#' @export
density.normal <- function(x, ...) {
    mu <- x$mu
    sd <- sqrt(x$var)
    function(t, log = FALSE, ...) {
        dnorm(x = t, mean = mu, sd = sd, log = log)
    }
}

#' Method for obtaining the cdf of an `normal` object.
#' @param x The object to obtain the cdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function \code{function(q, lower.tail = TRUE, log.p = FALSE, ...)}
#'   that computes the cdf (or log-cdf) of the normal distribution at \code{q}.
#' @examples
#' x <- normal(0, 1)
#' F <- cdf(x)
#' F(0)
#' F(1.96)
#' @importFrom stats pnorm
#' @export
cdf.normal <- function(x, ...) {
    mu <- x$mu
    sd <- sqrt(x$var)
    function(q, lower.tail = TRUE, log.p = FALSE, ...) {
        pnorm(q = q, mean = mu, sd = sd, lower.tail = lower.tail, log.p = log.p)
    }
}

#' Method for obtaining the support of a `normal` object, where the support
#' is defined as values that have non-zero probability density.
#' 
#' @param x The `normal` object to obtain the support of
#' @return A support-type object (see `support.R`), in this case an
#'         `interval` object for each component.
#' @examples
#' x <- normal(0, 1)
#' sup(x)
#' @export
sup.normal <- function(x) {
    interval$new()
}

#' Method for obtaining the dimension of a `normal` object.
#' @param x The `normal` object to obtain the dimension of
#' @return The dimension of the `normal` object
#' @examples
#' dim(normal(0, 1))
#' @export
dim.normal <- function(x) {
    1
}

#' Method for obtaining the inverse cdf of an `normal` object.
#' @param x The object to obtain the inverse cdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function \code{function(p, lower.tail = TRUE, log.p = FALSE, ...)}
#'   that computes the inverse cdf of the normal distribution.
#' @examples
#' x <- normal(0, 1)
#' q <- inv_cdf(x)
#' q(0.5)
#' q(0.975)
#' @importFrom stats qnorm
#' @export
inv_cdf.normal <- function(x, ...) {
    mu <- x$mu
    sd <- sqrt(x$var)
    function(p, lower.tail = TRUE, log.p = FALSE, ...) {
        qnorm(p, mean = mu, sd = sd, lower.tail = lower.tail, log.p = log.p)
    }
}


