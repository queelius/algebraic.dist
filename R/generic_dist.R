#' Generic method for obtaining the hazard function of an object.
#'
#' @param x The object to obtain the hazard function of.
#' @param ... Additional arguments to pass.
#' @export
hazard <- function(x, ...) {
    UseMethod("hazard", x)
}

#' Generic method for obtaining the pdf function of an object.
#'
#' @param x The object to obtain the pdf of.
#' @param ... Additional arguments to pass.
#' @export
pdf <- function(x, ...) {
    UseMethod("pdf", x)
}

#' Generic method for obtaining the cdf of an object.
#'
#' @param x The object to obtain the cdf of.
#' @param ... Additional arguments to pass.
#' @export
cdf <- function(x, ...) {
    UseMethod("cdf", x)
}

#' Generic method for obtaining the survival function of an object.
#'
#' @param x The object to obtain the survival function of.
#' @param ... Additional arguments to pass.
#' @export
surv <- function(x, ...) {
    UseMethod("surv", x)
}

#' Generic method for obtaining the quantile (inverse cdf) of an object.
#'
#' @param x The object to obtain the quantile of.
#' @param ... Additional arguments to pass.
#' @export
quantile <- function(x, ...) {
    UseMethod("quantile", x)
}


#' Generic method for obtaining the sampler of an object
#'
#' @param x The object to obtain the sampler of.
#' @param ... Additional arguments to pass.
#' @export
sampler <- function(x, ...) {
    UseMethod("sampler", x)
}

#' Generic method for obtaining the parameters of an object.
#' 
#' @param x The object to obtain the parameters of.
#' @param ... Additional arguments to pass.
#' @export
params <- function(x, ...) {
    UseMethod("params", x)
}

#' Generic method for obtaining the expectation of `f` with respect to
#' `x`.
#'
#' @param x The disrtibution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass.
#' @export
expectation <- function(x, g, ...) {
    UseMethod("expectation", x)
}

#' Generic method for obtaining the marginal distribution of a distribution
#' object `x` over components `indices`.
#' @param x The distribution object.
#' @param indices The indices of the marginal distribution to obtain.
#' @export
marginal <- function(x, indices) {
    UseMethod("marginal", x)
}


#' Generic method for obtaining the conditional distribution of a distribution
#' object `x` given condition `P`.
#' @param x The empirical distribution object.
#' @param P the predicate to condition on.
#' @export
conditional <- function(x, P) {
    UseMethod("conditional", x)
}


#' Generic method for applying a map `f` to distribution object `x`.
#' @param x The distribution object.
#' @param g The function to apply.
#' @param ... Additional arguments to pass.
#' @export
rmap <- function(x, g, ...) {
    UseMethod("rmap", x)
}


#' Generic method for retrieving the support of an object `x`.
#' 
#' The returned value should have the following operations:
#'  - `min`: a vector, the minimum value of the support for each component.
#'  - `max`: a vector, the maximum value of the support for each component.
#'  - `call`: a predicate function, which returns TRUE if the value is in
#'    the support, and FALSE otherwise.
#'  - `sample`: a function, which returns a sample from the support. Note that
#'    the returned value is not guaranateed to be in the support. You may need
#'    to call `call` to check.
#' @param x The object to obtain the support of.
#' @param ... Additional arguments to pass.
#' @return A support object for `x`.
#' @export
sup <- function(x, ...) {
    UseMethod("sup", x)
}