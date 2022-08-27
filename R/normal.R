#' Construct (multivariate or univariate) normal distribution object.
#'
#' @param mu mean
#' @param sigma variance-covariance matrix
#' @export
make_normal <- function(
    mu=c(0),
    sigma=diag(length(mu)))
{
    structure(.Data=c(mu=mu,sigma=sigma),class=c("normal","dist"))
}


print.normal <- function(x,...)
{
    if (x["mu"] == 0 && x["sigma"]==1)
    {
        cat("standard normal distribution\n")
    }
    else
    {
        cat("normal distribution\n")
        cat("parameters:\n")
        print.default(unclass(x),...)
    }
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
