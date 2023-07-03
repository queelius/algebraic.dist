#' Method for obtaining the expectation of `f` with respect to a
#' `univariate_dist` object `x`.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass into `pdf` and `sup`.
#' @param control An (optional) list of control parameters for `integrate`.
#' @importFrom stats integrate
#' @export
expectation.univariate_dist <- function(x, ..., control = list()) {
    # we assume the support is a contiguous interval that has operations
    # for retrieving the lower and upper bounds.
    S <- sup(x, ...)
    pdf.x <- pdf(x, ...)

    defaults <- list(
        lower = infimum(S),
        upper = supremum(S),
        subdivisions = 00L,
        rel.tol = 1e-20,
        abs.tol = 1e-20,
        stop.on.error = TRUE
    )
    control <- modifyList(defaults, control)
  
    function(
        g, ...,
        lower = control$lower,
        upper = control$upper,
        subdivisions = control$subdivisions,
        rel.tol = control$rel.tol,
        abs.tol = control$abs.tol,
        stop.on.error = control$stop.on.error) {

        r <- function(t, ...) g(t, ...) * pdf.x(t)
        integrate(
            f = r,
            lower = infimum(S),
            upper = supremum(S),
            ...,
            subdivisions = subdivisions,
            rel.tol = rel.tol,
            abs.tol = abs.tol,
            stop.on.error = stop.on.error)
    }
}


#' Method for obtaining the mean of `univariate_dist` object `x`.
#'
#' @param x The distribution object.
#' @param ... Additional arguments to pass.
#' @export
mean.univariate_dist <- function(x, ...) {
    expectation(x, function(t) t, ...)$value
}

#' Method for obtaining the variance of `univariate_dist` object `x`.
#'
#' @param x The distribution object.
#' @param ... Additional arguments to pass.
#' @export
var.univariate_dist <- function(x,...) {
    mu <- mean(x, ...)
    expectation(x, function(t) (t - mu)^2, ...)$value
}

#' Method for obtaining the standard deviation of `univariate_dist` object `x`.
#'
#' @param x The distribution object.
#' @param ... Additional arguments to pass.
#' @export
sd.univariate_dist <- function(x, ...) {
    sqrt(var(x, ...))
}

#' Method for obtaining the inverse cdf for a `univariate_dist` object.
#'
#' We use Newton's method to solve for the root of `cdf(x)(t) - p`.
#'
#' @param x The object to obtain the inverse cdf of.
#' @param ... Additional arguments to pass.
#' @export
inv_cdf.normal <- function(x, t0 = NULL, eps = 1e-3, ...) {
  # we assume that sup(x) is a contiguous interval
  if (is.null(t0)) {
      t0 <- sampler(x)(1)
  }
  stopifnot(has(Sx, t0))
  function(p, par, ...) {
    stopifnot(p >= 0 && p <= 1)

    par <- params(x, par)
    Sx <- sup(x, par, ...)
    Fx <- cdf(x, par, ...)
    fx <- pdf(x, par, ...)

    t1 <- NULL
    repeat {
      alpha <- 1
      repeat
      {
        t1 <- t0 - alpha * (Fx(t0, ...) - p) / fx(t0, ...)
        print(t1)
        if (has(Sx, t1)) {
          break
        }
        alpha <- alpha / 2
      }
      if (abs(t1 - t0) < eps) {
        break
      }
      t0 <- t1
    }
    t1
  }
}
