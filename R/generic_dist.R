#' Generic method for obtaining the hazard function of
#' an object.
#'
#' @param x The object to obtain the hazard function of
#' @param ... Additional arguments to pass
#' @export
hazard <- function(x, ...)
{
    UseMethod("hazard",x)
}

#' Generic method for obtaining the pdf function of
#' an object.
#'
#' @param x The object to obtain the hazard function of
#' @param ... Additional arguments to pass
#' @export
pdf <- function(x, ...)
{
    UseMethod("pdf",x)
}

#' Generic method for obtaining the survival function of
#' an object.
#'
#' @param x The object to obtain the survival function of
#' @param ... Additional arguments to pass
#' @export
surv <- function(x, ...)
{
    UseMethod("surv",x)
}

#' Generic method for obtaining the survival function of
#' a \code{dist} (distribution) object.
#'
#' @param x The object to obtain the cdf of
#' @param ... Additional arguments to pass
#' @export
cdf <- function(x, ...)
{
    UseMethod("cdf",x)
}
