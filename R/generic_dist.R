#' Generic method for obtaining the hazard function of an object.
#'
#' @param x The object to obtain the hazard function of.
#' @param ... Additional arguments to pass.
#' @export
hazard <- function(x, ...)
{
    UseMethod("hazard",x)
}

#' Generic method for obtaining the pdf function of an object.
#'
#' @param x The object to obtain the hazard function of.
#' @param ... Additional arguments to pass.
#' @export
pdf <- function(x, ...)
{
    UseMethod("pdf",x)
}

#' Generic method for obtaining the survival function of an object.
#'
#' @param x The object to obtain the survival function of.
#' @param ... Additional arguments to pass.
#' @export
surv <- function(x, ...)
{
    UseMethod("surv",x)
}

#' Generic method for obtaining the survival function of an object.
#'
#' @param x The object to obtain the cdf of.
#' @param ... Additional arguments to pass.
#' @export
cdf <- function(x, ...)
{
    UseMethod("cdf",x)
}

#' Generic method for obtaining the sampler of an object
#'
#' @param x The object to obtain the cdf of.
#' @param ... Additional arguments to pass.
#' @export
sampler <- function(x, ...)
{
    UseMethod("sampler",x)
}

#' Generic method for obtaining the inverse cdf of an object.
#'
#' @param x The object to obtain the inverse cdf of.
#' @param ... Additional arguments to pass.
#' @export
inv_cdf <- function(x,...)
{
    UseMethod("inv_cdf",x)
}

#' Generic method for obtaining the expectation of \code{f} with respect to
#' \code{x}.
#'
#' @param x The disrtibution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass.
#' @export
expectation <- function(x,g,...)
{
    UseMethod("expectation",x)
}
