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
#' @return The expected value (numeric scalar), or the full
#'   \code{integrate()} result if \code{compute_stats = TRUE}.
#' @examples
#' x <- normal(3, 4)
#' # E[X] for Normal(3, 4) is 3
#' expectation(x, function(t) t)
#'
#' # E[X^2] for Exp(1) is 2
#' expectation(exponential(1), function(t) t^2)
#' @importFrom stats integrate density
#' @importFrom utils modifyList
#' @export
expectation.univariate_dist <- function(x, g, ..., control = list()) {

    if (!("continuous_dist" %in% class(x))) {
        return(expectation.dist(x, g, ..., control))
    }

    S <- sup(x)
    f <- density(x)
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
#' @return Numeric scalar; the mean of the distribution.
#' @examples
#' mean(normal(5, 2))    # 5
#' mean(exponential(2))  # 0.5
#' @export
mean.univariate_dist <- function(x, ...) {
    expectation(x, function(t) t, ...)
}

#' Method for obtaining the variance of `univariate_dist` object.
#'
#' @param object The distribution object.
#' @param ... Additional arguments to pass into `expectation`.
#' @return Numeric scalar; the variance of the distribution.
#' @examples
#' vcov(normal(0, 4))    # 4
#' vcov(exponential(2))  # 0.25
#' @export
vcov.univariate_dist <- function(object, ...) {
    mu <- mean(object)
    expectation(object, function(t) (t - mu)^2, ...)
}
