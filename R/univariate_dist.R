#' \code{univariate_dist} models the concept of a continuous univariate
#' distribution.
#'
#' For an object that inherits \code{univariate_dist}, it will have default
#' implementations for everything *except* the \code{cdf} and the \code{sup}
#' (support) functions.
#'
#' We chose the cdf since it is easy to build each of the following from it:
#'     pdf, surv, inv_cdf, hazard, and sampler
#' The survival function (surv) would have also been an appropriate choice.
#'
#' Ideally, a subclass of a univariate distribution object will provide custom
#' implementations that are more efficient or accurate, but it is not required
#' and the defaults work well enough in most cases.
#' 
#' However, if the support is not contiguous, then \code{expectation} and
#' \code{inv_cdf} will need to be overridden from the default.
#' @export
NULL

univariate_dist <- function(cdf,sup=function(x) c(-Inf,Inf))
{
  structure(list(
    cdf=cdf,
    sup=sup),
    class="univariate_dist")
}

#' Method for obtaining the hazard function of a \code{univariate_dist} object.
#'
#' @param x The object to obtain the hazard function of
#' @param ... Additional arguments to pass
#' @export
hazard.univariate_dist <- function(x, ...)
{
  fx <- pdf(x,...)
  sx <- surv(x,...)
  function(t) fx(t)/sx(t)
}

#' Method for obtaining the pdf of a \code{univariate_dist} object.
#'
#' @param x The object to obtain the pdf of.
#' @param ... Additional arguments to pass.
#' @importFrom numDeriv grad
#' @export
pdf.univariate_dist <- function(x,...)
{
  Fx <- cdf(x,...)
  function(t) grad(Fx,t)
}

#' Method for obtaining the survival function a \code{univariate_dist} object.
#'
#' @param x The object to obtain the survival function of.
#' @param ... Additional arguments to pass.
#' @export
surv.univariate_dist <- function(x,...)
{
  Fx <- cdf(x,...)
  function(t) 1-Fx(t)
}

#' Method for obtaining the sampler for a \code{univariate_dist} object.
#'
#' We use the inverse cdf (\code{inv_cdf(x)}) and apply it to a uniform random
#' variable between 0 and 1.
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
#' We use Newton's method to solve for the root of \code{cdf(x)(t) - p}.
#'
#' @param x The object to obtain the inverse cdf of.
#' @param ... Additional arguments to pass.
#' @export
inv_cdf.univariate_dist <- function(x,t0=NULL,eps=1e-3,...)
{
  # we assume that sup(x) is a contiguous interval
  if (is.null(t0))
      t0 <- sampler(x)(1)
  sx <- sup(x)
  stopifnot(t0 %in% sx)
  Fx <- cdf(x,...)
  fx <- pdf(x,...)
  function(p)
  {
    stopifnot(p >= 0 && p <= 1)
    t1 <- NULL
    repeat
    {
      alpha <- 1
      repeat
      {
        t1 <- t0 - alpha * (Fx(t0)-p)/fx(t0)
        if (t1 %in% sx)
          break
        alpha <- alpha / 2
      }
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
    # we assume the support is a contiguous interval that has operations
    # for retrieving the lower and upper bounds.
    integrate(f=g,lower=lower_bound(sup(x)),upper=upper_bound(sup(x)),...)
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

#' Method for obtaining the standard deviation of \code{univariate_dist} object \code{x}.
#'
#' @param x The distribution object.
#' @param ... Additional arguments to pass.
#' @export
sd.univariate_dist <- function(x,...)
{
    sqrt(variance(x,...))
}

