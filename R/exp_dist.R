#' Construct exponential distribution object.
#'
#' @param rate failure rate
#' @export
exp_dist <- function(rate)
{
  structure(rate,
            class=unique(c("exp_dist","univariate_dist","dist",class(rate))))
}

#' Function to determine whether an object \code{x} is an \code{exp_dist} object.
#' @export
is_exp_dist <- function(x)
{
  inherits(x,"exp_dist")
}

#' @export
print.exp_dist <- function(x,...)
{
  cat("Exponential distribution with failure rate",x,"\n")
}

#' Method for obtaining the variance of a \code{exp_dist} object.
#'
#' @param object The \code{exp_dist} object to obtain the variance of
#' @param ... Additional arguments to pass
#' @importFrom stats vcov
#' @export
vcov.exp_dist <- function(object,...)
{
    1/unclass(object)^2
}

#' Method to obtain the hazard function of
#' an \code{exp_dist} object.
#'
#' @param x The \code{exp_dist} object to obtain the hazard function of
#' @param ... Additional arguments to pass
#' @export
hazard.exp_dist <- function(x,...)
{
    rate <- unclass(x)
    function(t) ifelse(t <= 0,0,rate)
}

#' Method to obtain the pdf of an \code{exp_dist} object.
#'
#' @param x The object to obtain the pdf of
#' @param logp Whether to compute the log of the density/prob
#' @export
pdf.exp_dist <- function(x,logp=F)
{
    rate <- unclass(x)
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
  rate <- unclass(x)
  function(n=1) as.matrix(stats::rexp(n,rate))
}


#' Method to obtain the minimum of a \code{exp_dist} object and another
#' distribution object. If the other distribution object is also an
#' \code{exp_dist} object, then the result is an \code{exp_dist} object
#' with the sum of the rates of the two \code{exp_dist} objects.
#' Otherwise, the result is the result of \code{min} applied to the
#' two distribution objects. We use \code{dispatch_dist} to dispatch
#' to the appropriate method for two objects of type \code{dist}.
#' 
#' @param x1 \code{exp_dist} object
#' @param x2 \code{dist} object
minimum.exp_dist <- function(x1,x2)
{
    if (is_exp_dist(x2))
      return(make_exp_dist(rate(x1)+rate(x2)))
    else
      return(dispatch(x1,x2,min))
}

