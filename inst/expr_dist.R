#' takes a function `f` and a distribution object `dist` and return a lazy
#' distribution that is only computed on demand.
#' @param f The function to apply.
#' @param dist The distribution object.
expr_dist <- function(f, dist) {
    structure(list(
        f = f,
        dist = dist),
        class = c(paste0(expr,"_dist"),"expr_dist"))
}

#' Function to determine whether an object `x` is an `expr_dist` object.
#' @export
#' @param x The object to test
is_expr_dist <- function(x) {
    inherits(x, "expr_dist")
}


#' Method for obtaining the sampler of an `expr_dist` object.
#' @param x The `expr_dist` object to obtain the sampler of.
#' @param n The number of samples to obtain.
#' @param ... Additional arguments to pass.
#' @export
sampler.expr_dist <- function(x, ...) {
    rx <- sampler(x$dist)
    function(n = 1) {
        x$f(rx(n, ...))
    }
}

#' Method for obtaining the expectation of an `expr_dist` object.
#' @param x The `expr_dist` object to obtain the expectation of.
#' @param ... Additional arguments to pass.
#' @export
expectation.expr_dist <- function(x, f, ...) {

    rx <- sampler(x$dist)
    function(n = 1) {
        mean(x$f(rx(n, ...)))
    }
}

#' Method for obtaining the variance of an `expr_dist` object.
#' @param x The `expr_dist` object to obtain the variance of.
#' @param ... Additional arguments to pass.
#' @export
#' @importFrom stats vcov
vcov.expr_dist <- function(x, ...) {
    rx <- sampler(x$dist)
    function(n = 1) {
        cov(x$f(rx(n, ...)))
    }
}

#' Method for obtaining the pdf of an `expr_dist` object.
#' @param x The `expr_dist` object to obtain the pdf of.
#' @param ... Additional arguments to pass.
#' @export
pdf.expr_dist <- function(x, ...) {

    # we use the sampler and estimate the kernel density
    # of the function
    rx <- sampler(x$dist)
    kde(x$f(rx(n, ...)))
}