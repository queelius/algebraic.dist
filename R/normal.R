#' Construct (multivariate or univariate) normal distribution object.
#'
#' @param mu mean
#' @param sigma variance-covariance matrix
#' @export
make_normal <- function(
    mu=c(0),
    sigma=diag(length(mu)))
{
    structure(list(
        mu=mu,
        sigma=sigma),
        class=c("normal","dist"))
}

#' Retrieve the variance-covariance matrix (or scalar)
#' of a \code{normal} object.

#' A variance-covariance matrix is a square matrix
#' giving the covariance between each pair of elements
#' of a given random vector. Importantly, its diagonal
#' contains variances of the elements
#'
#' @param object The \code{normal} object to retrieve
#'               the variance-covariance matrix from
#' @export
vcov.normal <- function(object,...)
{
    object$sigma
}

#' Retrieve the fisher information of a \code{normal} object for
#' parameters \code{mu} and \code{sigma}.

#' Fisher information is a way of measuring the amount of
#' information that an observable random variable `X`
#' carries about an unknown parameter \code{theta}
#' upon which the probability of `X` depends.
#'
#' The inverse of the Fisher information matrix
#' is the variance-covariance of the MLE for
#' \code{theta}.
#'
#' @param object The \code{normal} object
#' @export
fisher_info.normal <- function(object,...)
{
    # incorrect, i think, revisit
    matrix(1/object$sigma,1/(2*object$sigma^2))
}

#' Method for obtaining the parameters of a \code{normal} object.
#'
#' @param x The object to obtain the parameters of
#' @importFrom algebraic.mle params
#' @export
params.normal <- function(x, ...)
{
    c("mu"=x$mu,"sigma"=x$sigma)
}

#' Method for sampling from a \code{normal} object.
#'
#' @param x The object to sample from
#' @importFrom mvtnorm rmvnorm
#' @importFrom algebraic.mle sampler
#' @export
sampler.normal <- function(x, ...)
{
    function(n=1)
    {
        mvtnorm::rmvnorm(n,x$mu,x$sigma)
    }
}
