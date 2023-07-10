#' Method for obtaining the expectation of `f` with respect to a
#' `univariate_dist` object `x`.
#' 
#' Assumes the support is a contiguous interval that has operations
#' for retrieving the lower and upper bounds.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass into `g`.
#' @param control An (optional) list of control parameters for `integrate` or
#'               `expectation_data` (if `x` is not continuous)
#' @importFrom stats integrate
#' @importFrom utils modifyList
#' @export
expectation.univariate_dist <- function(x, g, ..., control = list()) {

    if (!("continuous_dist" %in% class(x))) {
        return(expectation.dist(x, g, ..., control))
    }

    S <- sup(x)
    f <- pdf(x)
    defaults <- formals(integrate)
    defaults$compute_stats <- FALSE
    defaults$rel.tol <- 1e-4
    defaults$abs.tol <- defaults$rel.tol
    control <- modifyList(defaults, control)

    res <- integrate(
        function(t) g(t, ...) * f(t),
        lower = infimum(S),
        upper = supremum(S),
        subdivisions = control$subdivision,
        rel.tol = control$rel.tol,
        abs.tol = control$abs.tol,
        stop.on.error = control$stop.on.error)

    if (control$compute_stats) {
        return(res)
    } else {
        return(res$value)
    }
}

#' Method for obtaining the mean of `univariate_dist` object `x`.
#'
#' @param x The distribution object.
#' @param ... Additional arguments to pass into `expectation`.
#' @export
mean.univariate_dist <- function(x, ...) {
    expectation(x, function(t) t, ...)
}

#' Method for obtaining the variance of `univariate_dist` object `x`.
#'
#' @param x The distribution object.
#' @export
vcov.univariate_dist <- function(x) {
    mu <- mean(x)
    expectation(x, function(t) (t - mu)^2)
}

#' Method for obtaining the standard deviation of `univariate_dist` object `x`.
#'
#' @param x The distribution object.
#' @export
sd.univariate_dist <- function(x) {
    sqrt(vcov(x))
}


