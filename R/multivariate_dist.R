#' \code{multivariate_dist} models the concept of a continuous multivariate
#' distribution.
#' @export
NULL

#' Method for obtaining the pdf of a \code{univariate_dist} object.
#'
#' @param x The object to obtain the pdf of.
#' @param ... Additional arguments to pass.
#' @importFrom numDeriv grad
#' @export
pdf.multivariate_dist <- function(x,...)
{
    cdfx <- cdf(x,...)
    function(t) grad(cdfx,t)
}

#' Method for obtaining the sampler for a \code{multivariate_dist} object.
#'
#' @param x The object to obtain the sampler of.
#' @param ... Additional arguments to pass.
#' @export
sampler.univariate_dist <- function(x,...)
{
    q <- inv_cdf(x,...)
    function(n=1) q(runif(n))
}

#' Method for obtaining the inverse cdf for a \code{univariate_dist} object.
#'
#' We use Newton's method to solve for the root of cdf(x)(t) - p.
#'
#' @param x The object to obtain the inverse cdf of.
#' @param ... Additional arguments to pass.
#' @export
inv_cdf.univariate_dist <- function(x,t0=NULL,eps=1e-3,...)
{
    if (is.null(t0))
        t0 <- sampler(x)(1)
    stopifnot(sup(x)(t0))
    cdfx <- cdf(x,...)
    pdfx <- pdf(x,...)
    function(p)
    {
        stopifnot(p >= 0 && p <= 1)
        t1 <- NULL
        repeat
        {
            t1 <- t0 - (cdfx(t0)-p)/pdfx(t0)
            if (abs(t1-t0) < eps)
                break
            t0 <- t1
        }
        t1
    }
}

#' Method for obtaining the expectation of \code{f} with respect to a
#' \code{univariate_dist} object \code{x}.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass.
#' @importFrom stats integrate
#' @export
expectation.univariate_dist <- function(x,g,...)
{
    # we assume the support is a contiguous interval, represented by
    # a vector where s[1] and s[2] are respectively the lower and upper bounds
    s <- sup(x)
    integrate(f=g,lower=s[1],upper=s[2],...)
}

#' Method for obtaining the mean of \code{univariate_dist} object \code{x}.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass.
#' @importFrom stats integrate
#' @export
mean.univariate_dist <- function(x,...)
{
    expectation(x,function(t) t,...)
}

#' Method for obtaining the variance of \code{univariate_dist} object \code{x}.
#'
#' @param x The distribution object.
#' @param ... Additional arguments to pass.
#' @export
variance.univariate_dist <- function(x,...)
{
    mu <- mean(x,...)
    expectation(x,function(t) (t-mu)^2,...)
}
