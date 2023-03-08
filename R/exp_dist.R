#' Construct exponential distribution object.
#'
#' @param rate failure rate
#' @export
exp_dist <- function(rate)
{
    structure(rate,
              class=unique(c("exp_dist","univariate_dist","dist",class(rate))))
}

#' @export
is_exp_dist <- function(x)
{
  inherits(x,"exp_dist")
}

#' @export
print.exp_dist <- function(x,...)
{
  cat("Exponential distribution with failure rate",unclass(x))
}

#' Method for obtaining the variance of a \code{exp_dist} object.
#'
#' @param object The \code{exp_dist} object to obtain the variance of
#' @param ... Additional arguments to pass
#' @importFrom stats vcov
#' @export
vcov.exp_dist <- function(object,...)
{
    1/params(x)^2
}

#' Method for obtaining the parameters of a \code{exp_dist} distribution object.
#'
#' @param x The \code{exp_dist} object to obtain the parameters of.
# #' @importFrom algebraic.mle params
#' @export
params.exp_dist <- function(x)
{
    unclass(x)
}

#' Method to obtain the hazard function of
#' an \code{exp_dist} object.
#'
#' @param x The \code{exp_dist} object to obtain the hazard function of
#' @param ... Additional arguments to pass
#' @export
hazard.exp_dist <- function(x,...)
{
    rate <- params(x)
    function(t) ifelse(t <= 0,0,rate)
}

#' Method to obtain the pdf of an \code{exp_dist} object.
#'
#' @param x The object to obtain the pdf of
#' @param logp Whether to compute the log of the density/prob
#' @export
pdf.exp_dist <- function(x,logp=F)
{
    rate <- params(x)
    ifelse(logp,
           function(t) ifelse(t<=0,-Inf,log(rate) - rate*t),
           function(t) ifelse(t<=0,0,rate*exp(-rate*t)))
}

#' Method to sample from an \code{exp_dist} object.
#'
#' @param x The \code{exp_dist} object to sample from.
#' @importFrom algebraic.mle sampler
#' @export
sampler.exp_dist <- function(x)
{
  rate <- params(x)
  function(n=1) as.matrix(stats::rexp(n,rate))
}

