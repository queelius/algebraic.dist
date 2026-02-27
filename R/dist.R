#' Function to determine whether an object `x` is a `dist` object.
#' @param x The object to test
#' @return Logical indicating whether `x` is a `dist` object.
#' @export
is_dist <- function(x) {
  inherits(x, "dist")
}

#' Expectation of a Function Applied to a `dist` Object
#'
#' Expectation operator applied to `x` of type `dist`
#' with respect to a function `g`. Optionally, constructs a confidence interval
#' for the expectation estimate using the Central Limit Theorem.
#'
#' @param x A `dist` object.
#' @param g Characteristic function of interest, defaults to identity. 
#' @param ... Additional arguments to pass to `g`.
#' @param control A list of control parameters:
#'  compute_stats - Logical, whether to compute CIs for the expectations, defaults
#'                  to FALSE
#'  n             - Integer, the number of samples to use for the MC estimate,
#'                  defaults to 10000L
#'  alpha         - Real, the significance level for the confidence interval,
#'                  defaults to 0.05
#' @return If `compute_stats` is FALSE, then the estimate of the expectation,
#'        otherwise a list with the following components:
#'   value - The estimate of the expectation
#'   ci    - The confidence intervals for each component of the expectation
#'   n     - The number of samples
#' @export
expectation.dist <- function(
    x,
    g = function(t) t,
    ...,
    control = list()) {

    defaults <- list(
      compute_stats = FALSE,
      n = 10000L,
      alpha = 0.05)
    control <- modifyList(defaults, control)

    stopifnot(is.numeric(control$n), control$n > 0,
      is.numeric(control$alpha), control$alpha > 0, control$alpha < 1)
  
    expectation_data(
      data = sampler(x)(control$n),
      g = g,
      ...,
      compute_stats = control$compute_stats,
      alpha = control$alpha)
}

#' Method for obtaining a summary of a `dist` object.
#' @param object The object to obtain the summary of
#' @param ... Additional arguments to pass
#' @param name The name of the distribution, defaults to the class of the object.
#' @param nobs The number of observations to report for the summary, if applicable.
#' @return A `summary_dist` object
#' @export
summary.dist <- function(object, ..., name = NULL, nobs = NULL) {
  summary_dist(
    name = if (is.null(name)) class(object)[1] else name,
    mean = mean(object),
    vcov = vcov(object),
    nobs = nobs)
}

#' Method for obtaining the condition distribution, `x | P(x)`, of
#' `dist` object `x`.
#' 
#' Falls back to MC: materializes `x` via `ensure_realized()` and
#' then conditions on the resulting empirical distribution.
#' 
#' @param x The distribution object.
#' @param P The predicate function to condition the distribution on
#' @param n The number of samples to generate for the MC estimate of the
#'         conditional distribution x | P. Defaults to 10000.
#' @param ... additional arguments to pass into `P`.
#' @export
conditional.dist <- function(x, P, n = 10000L, ...) {
  conditional(ensure_realized(x, n = n), P, ...)
}

#' Method for obtaining g(x)) where x is a `dist` object.
#' 
#' Falls back to MC: materializes `x` via `ensure_realized()` and
#' then applies `rmap` with `g` to the resulting empirical distribution.
#' @param x The distribution object.
#' @param g The function to apply to the distribution.
#' @param n The number of samples to generate for the MC estimate of the
#'          conditional distribution x | P. Defaults to 10000.
#' @param ... additional arguments to pass into `g`.
#' @export
rmap.dist <- function(x, g, n = 10000L, ...) {
  rmap(ensure_realized(x, n = n), g, ...)
}


#' Print method for `summary_dist` objects.
#' @param x The object to print
#' @param ... Additional arguments
#' @export
print.summary_dist <- function(x, ...) {
  cat(x$name, "\n")
  cat("Mean:\n")
  print(x$mean)
  cat("Covariance:\n")
  print(x$vcov)
  if (!is.null(x$nobs)) {
    cat("Number of observations:", x$nobs, "\n")
  }
}

#' Method for constructing a `summary_dist` object.
#' @param name The name of the distribution
#' @param mean The mean of the distribution
#' @param vcov The variance of the distribution
#' @param nobs The number of observations used to construct the distribution,
#'             if applicable.
#' @return A `summary_dist` object
#' @export
summary_dist <- function(name, mean, vcov, nobs = NULL) {
  structure(list(
    name = name,
    mean = mean,
    vcov = vcov,
    nobs = nobs),
    class = "summary_dist")
}





#' Sampler for non-dist objects (degenerate distributions).
#'
#' @param x The object to sample from
#' @param ... Additional arguments to pass
#' @return A function that takes n and returns n copies of x
#' @export
sampler.default <- function(x, ...) {
  function(n) {
    rep(x, n)
  }
}

#' Variance-covariance for non-dist objects (degenerate distributions).
#'
#' @param object The object (returns 0 for constants)
#' @param ... Additional arguments to pass (not used)
#' @return 0 (degenerate distributions have no variance)
#' @export
vcov.default <- function(object, ...) {
  0
}

