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
#' @examples
#' # Bivariate normal with identity covariance
#' X <- mvn(mu = c(0, 0))
#' mean(X)
#' vcov(X)
#'
#' # 1D case returns a normal object
#' is_normal(mvn(mu = 1, sigma = matrix(4)))
#' @export
mvn <- function(
    mu,
    sigma = diag(length(mu))) {

    if (!is.numeric(mu))
        stop("'mu' must be a numeric vector, got: ", deparse(mu))
    if (!is.matrix(sigma) || !is.numeric(sigma))
        stop("'sigma' must be a numeric matrix")
    if (nrow(sigma) != ncol(sigma))
        stop("'sigma' must be square, got ", nrow(sigma), "x", ncol(sigma))
    if (nrow(sigma) != length(mu))
        stop("dimensions of 'sigma' (", nrow(sigma), ") must match length of 'mu' (", length(mu), ")")

    if (length(mu) == 1) {
        return(normal(mu, sigma))
    }
    structure(list(mu = mu, sigma = sigma),
        class = c("mvn", "multivariate_dist",
                  "continuous_dist", "dist"))
}

#' Retrieve the variance-covariance matrix of an `mvn` object.
#' @param object The `mvn` object to retrieve the variance-covariance matrix of
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `mvn` object
#' @examples
#' X <- mvn(c(0, 0), diag(2))
#' vcov(X)
#' @export
vcov.mvn <- function(object, ...) {
    object$sigma
}

#' Retrieve the mean of a `mvn` object.
#' @param x The `mvn` object to retrieve the mean from
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `mvn` object
#' @examples
#' X <- mvn(c(1, 2, 3))
#' mean(X)
#' @export
mean.mvn <- function(x, ...) {
    x$mu
}

#' Method for obtaining the parameters of a `mvn` object.
#'
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @examples
#' X <- mvn(c(0, 0), diag(2))
#' params(X)
#' @export
params.mvn <- function(x) {
    c("mu" = x$mu, "sigma" = x$sigma)
}

#' Function to determine whether an object `x` is an `mvn` object.
#' @param x The object to test
#' @return Logical; \code{TRUE} if \code{x} is an \code{mvn}.
#' @examples
#' is_mvn(mvn(c(0, 0)))
#' is_mvn(normal(0, 1))
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
#' @examples
#' X <- mvn(c(0, 0), diag(2))
#' s <- sampler(X)
#' set.seed(42)
#' s(3)
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
#' @examples
#' X <- mvn(c(0, 0), diag(2))
#' f <- density(X)
#' f(c(0, 0))
#' f(c(1, 1))
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
#' @examples
#' X <- mvn(c(0, 0))
#' sup(X)
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
#' @return A \code{normal} (for a single index) or \code{mvn} marginal distribution.
#' @examples
#' X <- mvn(c(1, 2, 3))
#' # Univariate marginal
#' marginal(X, 1)
#' # Bivariate marginal
#' marginal(X, c(1, 3))
#' @export
marginal.mvn <- function(x, indices) {
    if (length(indices) == 0)
        stop("indices must be non-empty")
    if (any(indices < 1) || any(indices > dim(x)))
        stop("indices must be in [1, dim(x)]")

    if (length(indices) == 1)
        return(normal(x$mu[indices], x$sigma[indices, indices]))

    mvn(x$mu[indices], x$sigma[indices, indices])
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
#' @return An \code{n} by \code{length(mu)} matrix of samples within the
#'   probability region.
#' @examples
#' \donttest{
#' set.seed(42)
#' pts <- sample_mvn_region(10, mu = c(0, 0), sigma = diag(2), p = 0.95)
#' dim(pts)
#' }
#' @importFrom stats qchisq mahalanobis
#' @export
sample_mvn_region <- function(n, mu, sigma, p = .95, ...) {

    stopifnot(p >= 0 && p <= 1, n > 0)

    k <- length(mu)
    crit <- qchisq(p, k)
    i <- 0L
    samp <- sampler(mvn(mu = mu, sigma = sigma))
    data <- matrix(nrow = n, ncol = k)
    while (i < n) {
        candidate <- samp(1)
        d <- mahalanobis(candidate, center = mu, cov = sigma, ...)
        if (d <= crit) {
            i <- i + 1L
            data[i, ] <- candidate
        }
    }
    data
}


#' Format method for `mvn` objects.
#' @param x The object to format
#' @param ... Additional arguments (not used)
#' @return A character string
#' @examples
#' format(mvn(c(0, 0)))
#' @export
format.mvn <- function(x, ...) {
    sprintf("Multivariate normal distribution (%d dimensions)", length(x$mu))
}

#' Method for printing an `mvn` object.
#' @param x The object to print
#' @param ... Additional arguments to pass to `print`
#' @return \code{x}, invisibly.
#' @examples
#' print(mvn(c(0, 0)))
#' @export
print.mvn <- function(x, ...) {
    cat(format(x), "\n")
    cat("  mu:\n")
    print(x$mu)
    cat("  sigma:\n")
    print(x$sigma)
    invisible(x)
}


#' Method for obtaining the dimension of an `mvn` object.
#' @param x The object to obtain the dimension of
#' @return The dimension of the `mvn` object
#' @examples
#' dim(mvn(c(0, 0, 0)))
#' @export
dim.mvn <- function(x) {
    length(x$mu)
}

#' Method for obtaining the CDF of a `mvn` object.
#' @param x The object to obtain the CDF of
#' @param ... Additional arguments to pass (not used)
#' @return A function computing the multivariate normal CDF.
#' @examples
#' X <- mvn(c(0, 0), diag(2))
#' F <- cdf(X)
#' F(c(0, 0))
#' @importFrom mvtnorm pmvnorm
#' @export
cdf.mvn <- function(x, ...) {
    function(q, mu = x$mu, sigma = x$sigma, lower.tail = TRUE,
             log.p = FALSE, ...) {
        p <- pmvnorm(mean = mu, sigma = sigma, upper = q, ...)
        if (attr(p, "error") != 0) {
            warning(sprintf("error (%s): %s", attr(p, "error"), attr(p, "msg")))
        }
        p <- as.numeric(p)
        if (log.p) log(p) else p
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
#' @return An \code{mvn} distribution fitted to the transformed samples.
#' @examples
#' \donttest{
#' X <- mvn(c(1, 2), diag(2))
#' set.seed(42)
#' Y <- rmap(X, function(x) x^2)
#' mean(Y)
#' }
#' @importFrom stats vcov
#' @export
rmap.mvn <- function(x, g, n = 10000L, ...) {
    D <- rmap(ensure_realized(x, n = n), g, ...)
    mvn(mu = mean(D), sigma = vcov(D))
}


#' Conditional distribution for multivariate normal.
#'
#' Supports two calling patterns:
#' \enumerate{
#'   \item \strong{Closed-form} (via \code{given_indices} and
#'     \code{given_values}): Uses the exact Schur complement formula.
#'     Returns a \code{normal} (1D result) or \code{mvn}.
#'   \item \strong{Predicate-based} (via \code{P}): Falls back to MC
#'     realization via \code{\link{ensure_realized}}.
#' }
#'
#' @param x An \code{mvn} object.
#' @param P Optional predicate function for MC fallback.
#' @param ... Additional arguments forwarded to the predicate \code{P}.
#' @param given_indices Integer vector of observed variable indices.
#' @param given_values Numeric vector of observed values (same length as
#'   \code{given_indices}).
#' @return A \code{normal}, \code{mvn}, or \code{empirical_dist} object.
#' @examples
#' # Closed-form conditioning: X2 | X1 = 1
#' sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' X <- mvn(c(0, 0), sigma)
#' X2_given <- conditional(X, given_indices = 1, given_values = 1)
#' mean(X2_given)
#' vcov(X2_given)
#'
#' # Predicate-based MC fallback (slower)
#' \donttest{
#' set.seed(42)
#' X2_mc <- conditional(X, P = function(x) x[1] > 0)
#' }
#' @export
conditional.mvn <- function(x, P = NULL, ...,
                            given_indices = NULL, given_values = NULL) {
    # Closed-form via Schur complement
    if (!is.null(given_indices) && !is.null(given_values)) {
        if (length(given_indices) != length(given_values))
            stop("'given_indices' and 'given_values' must have the same length")
        if (any(given_indices < 1) || any(given_indices > dim(x)))
            stop("'given_indices' must be in [1, dim(x)]")

        all_idx <- seq_len(dim(x))
        free_idx <- setdiff(all_idx, given_indices)
        if (length(free_idx) == 0)
            stop("cannot condition on all variables")

        mu1 <- x$mu[free_idx]
        mu2 <- x$mu[given_indices]
        sig11 <- x$sigma[free_idx, free_idx, drop = FALSE]
        sig12 <- x$sigma[free_idx, given_indices, drop = FALSE]
        sig22 <- x$sigma[given_indices, given_indices, drop = FALSE]

        sig22_inv <- solve(sig22)
        mu_cond <- as.numeric(mu1 + sig12 %*% sig22_inv %*% (given_values - mu2))
        sig_cond <- sig11 - sig12 %*% sig22_inv %*% t(sig12)

        if (length(free_idx) == 1) {
            return(normal(mu = mu_cond, var = as.numeric(sig_cond)))
        }
        return(mvn(mu = mu_cond, sigma = sig_cond))
    }

    # Predicate-based MC fallback
    if (!is.null(P)) {
        return(conditional(ensure_realized(x), P, ...))
    }

    stop("must provide either 'P' or both 'given_indices' and 'given_values'")
}


#' Affine transformation of a normal or multivariate normal distribution.
#'
#' Computes the distribution of \eqn{AX + b} where \eqn{X \sim MVN(\mu, \Sigma)}.
#' The result is \eqn{MVN(A\mu + b, A \Sigma A^T)}.
#'
#' For a univariate \code{normal}, scalars \code{A} and \code{b} are promoted
#' to 1x1 matrices and scalar internally. Returns a \code{normal} if the
#' result is 1-dimensional.
#'
#' @param x A \code{normal} or \code{mvn} distribution object.
#' @param A A numeric matrix (or scalar for univariate).
#' @param b An optional numeric vector (or scalar) for the offset. Default is
#'   a zero vector.
#' @return A \code{normal} or \code{mvn} distribution.
#' @examples
#' X <- mvn(c(0, 0), diag(2))
#' # Project to first component via 1x2 matrix
#' Y <- affine_transform(X, A = matrix(c(1, 0), 1, 2), b = 5)
#' mean(Y)
#'
#' # Scale a univariate normal
#' Z <- affine_transform(normal(0, 1), A = 3, b = 2)
#' mean(Z)
#' vcov(Z)
#' @export
affine_transform <- function(x, A, b = NULL) {
    if (!inherits(x, "dist"))
        stop("'x' must be a 'dist' object")

    if (is_normal(x)) {
        mu <- mean(x)
        sigma <- matrix(vcov(x), 1, 1)
    } else if (is_mvn(x)) {
        mu <- x$mu
        sigma <- x$sigma
    } else {
        stop("'x' must be a 'normal' or 'mvn' distribution")
    }

    if (!is.matrix(A))
        A <- matrix(A, nrow = 1, ncol = length(mu))
    if (ncol(A) != length(mu))
        stop("ncol(A) must equal dim(x), got ", ncol(A), " vs ", length(mu))

    if (is.null(b))
        b <- rep(0, nrow(A))
    if (length(b) != nrow(A))
        stop("length(b) must equal nrow(A), got ", length(b), " vs ", nrow(A))

    mu_new <- as.numeric(A %*% mu + b)
    sigma_new <- A %*% sigma %*% t(A)

    if (length(mu_new) == 1) {
        normal(mu = mu_new, var = as.numeric(sigma_new))
    } else {
        mvn(mu = mu_new, sigma = sigma_new)
    }
}
