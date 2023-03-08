#' Construct exponential distribution object.
#'
#' @param rate failure rate
#' @export
make_exp_dist <- function(rate)
{
    structure(list(theta=rate),
              class=c("exp_dist","univariate_dist","dist"))
}

#' Method for obtaining the variance of a \code{exp_dist} object.
#'
#' @param object The \code{exp_dist} object to obtain the variance of
#' @param ... Additional arguments to pass
#' @importFrom stats vcov
#' @export
vcov.exp_dist <- function(object, ...)
{
    1/params.exp_dist(x)^2
}

#' Method for obtaining the parameters of a \code{exp_dist} distribution object.
#'
#' @param x The \code{exp_dist} object to obtain the parameters of.
#' @param ... Additional arguments to pass.
# #' @importFrom algebraic.mle params
#' @export
params.exp_dist <- function(x,...)
{
    x$theta
}

#' Method to obtain the hazard function of
#' an \code{exp_dist} object.
#'
#' @param x The \code{exp_dist} object to obtain the hazard function of
#' @param ... Additional arguments to pass
#' @export
hazard.exp_dist <- function(x,...)
{
    theta <- params.exp_dist(x)
    function(t) ifelse(t <= 0,0,theta)
}

#' Method to obtain the pdf of an \code{exp_dist} object.
#'
#' @param x The object to obtain the pdf of
#' @export
pdf.exp_dist <- function(x,...)
{
    theta <- params(x)
    function(t) ifelse(t <= 0,0,theta*exp(-theta*t))
}

#' Method to sample from an \code{exp_dist} object.
#'
#' @param x The \code{exp_dist} object to sample from.
#' @importFrom algebraic.mle sampler
#' @export
sampler.exp_dist <- function(x,...)
{
    theta <- params(x)
    function(n=1) stats::rexp(n,theta,...)
}




min.exp_dist <- function(x1,x2)
{
    if (is_exp_dist(x2))
    make_exp_dist()
}
