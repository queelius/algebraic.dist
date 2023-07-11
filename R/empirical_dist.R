#' Construct empirical distribution object.
#' @param data data to construct empirical distribution from. if matrix or
#'             data frame, each row is a joint observation, if a vector, each
#'             element is an observation. whatever data is, it must be
#'             convertible to a tibble.
#' @export
empirical_dist <- function(data) {

  if (length(data) == 0) {
    stop("data must have at least one observation")
  }
  if (is.matrix(data)) {
    cls_names <- c("empirical_dist", "multivariate_dist",
                   "discrete_dist", "dist")
  } else {
    data = as.matrix(data)
    cls_names <- c("empirical_dist", "univariate_dist",
                   "discrete_dist", "dist")
  }
  structure(list(data = data), class = cls_names)
}

#' Function to determine whether an object `x` is an `empirical_dist` object.
#' @param x The object to test
#' @export
is_empirical_dist <- function(x) {
  inherits(x, "empirical_dist")
}

#' Method for obtaining the dimension of a `empirical_dist` object.
#' @param x The object to obtain the dimension of.
#' @export
dim.empirical_dist <- function(x) {
  ncol(x$data)
}

#' Method for obtaining the pdf of a `empirical_dist` object.
#'
#' @param x The object to obtain the pdf of.
#' @param ... Additional arguments to pass into the pdf function.
#' @note sort tibble lexographically and do a binary search to find upper
#'       and lower bound in `log(nobs(x))` time.
#' @export
density.empirical_dist <- function(x, ...) {
  n <- nobs(x)
  xobs <- obs(x)
  p <- dim(x)
  function(t, log = FALSE) {
    stopifnot(is.numeric(t))

    # compute number of elements (or rows) in obs(x) equal to t
    if (is.matrix(t)) {
      stopifnot(ncol(t) == p)
      counts <- apply(t, 1, function(t_ob) { sum(apply(xobs, 1, function(xob) {
        all(xob == t_ob)
      })) })
    } else {
      stopifnot(length(t) == p)
      counts <- sum(apply(xobs, 1, function(xob) {
        all(xob == t)
      }))
    }
    if (log) {
      return(log(counts) - log(n))
    } else {
      return(counts / n)
    }
  }
}

#' Method for obtaining the sampler for a `empirical_dist` object.
#'
#' @param x The object to obtain the sampler of.
#' @param ... Additional arguments to pass (not used).
#' @importFrom stats nobs
#' @export
sampler.empirical_dist <- function(x, ...) {
  R <- nobs(x, ...)
  stopifnot(R > 0)

  function (n = 1, ...) {
    x$data[sample(1:R, size = n, replace = TRUE, ...), ]
  }
}

#' Method for obtaining the expectation of `empirical_dist` object `x`
#' under function `g`.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass itno function `g`.
#' @param control a list of control parameters:
#'  compute_stats - Whether to compute CIs for the expectations, defaults
#'                  to FALSE
#'  n             - The number of samples to use for the MC estimate,
#'                  defaults to 10000
#'  alpha         - The significance level for the confidence interval,
#'                  defaults to 0.05
#' @return If `compute_stats` is FALSE, then the estimate of the expectation,
#'        otherwise a list with the following components:
#'   value - The estimate of the expectation
#'   ci    - The confidence intervals for each component of the expectation
#'   n     - The number of samples
#' @export
expectation.empirical_dist <- function(
  x,
  g = function(t) t,
  ...,
  control = list()) {

    defaults <- list(
      compute_stats = FALSE,
      alpha = 0.05)
    control <- modifyList(defaults, control)

    stopifnot(is.numeric(control$alpha), control$alpha > 0, control$alpha < 1)
  
    expectation_data(
      data = obs(x),
      g = g,
      ...,
      compute_stats = control$compute_stats,
      alpha = control$alpha)
  }

#' Method for obtaining the mean of `empirical_dist` object `x`.
#'
#' @param x The distribution object.
#' @param ... Additional arguments to pass (not used).
#' @export
mean.empirical_dist <- function(x, ...) {
  colMeans(x$data)
}

#' Method for obtaining the variance of `empirical_dist` object `x`.
#'
#' @param object The empirical distribution object.
#' @param ... Additional arguments to pass (not used).
#' @importFrom stats cov
#' @export
vcov.empirical_dist <- function(object, ...) {
  cov(object$data)
}

#' Method for obtaining the marginal distribution of `empirical_dist` object
#' `x`.
#' @param x The empirical distribution object.
#' @param indices The indices of the marginal distribution to obtain.
#' @export
marginal.empirical_dist <- function(x, indices) {
  if (any(indices < 1) || any(indices > dim(x))) {
    stop("Invalid indices")
  }
  empirical_dist(x$data[, indices])
}

#' Method for obtaining the condition distribution, `x | P(x)`, of
#' `empirical_dist` object `x`.
#' 
#' In other words, we condition the data on the predicate function.
#' In order to do so, we simply remove all rows from the data that do not
#' satisfy the predicate P. For instance, if we have a 2-dimensional
#' distribution, and we want to condition on the first dimension being
#' greater than the second dimension, we would do the following:
#' 
#' ```r
#' x_cond <- conditional(x, function(d) d[1] > d[2])
#' ```
#' 
#' This would return a new empirical distribution object with the same
#' dimensionality as `x`, but with all rows where the first dimension is
#' less than or equal to the second dimension removed. 
#' 
#' @param x The empirical distribution object.
#' @param P The predicate function to condition the data on.
#' @param ... additional arguments to pass into `P`.
#' @export
conditional.empirical_dist <- function(x, P, ...) {
  empirical_dist(x$data[apply(X = x$data, MARGIN = 1, FUN = P, ...), ])
}


#' Method for obtaining the empirical distribution of a function of the
#' observations of `empirical_dist` object `x`.
#' 
#' @param x The empirical distribution object.
#' @param g The function to apply to each observation.
#' @param ... Additional arguments to pass into function `g`.
#' @export
rmap.empirical_dist <- function(x, g, ...) {
  g_data <- apply(x$data, 1, g, ...)
  if (is.matrix(g_data)) {
    g_data <- t(g_data)
  }
  empirical_dist(g_data)
}

#' Method for obtaining the cdf of `empirical_dist` object `x`.
#'
#' If `x` is a multivariate empirical distribution, this function will
#' throw an error. It's only defined for univariate empirical distributions.
#' @param x The empirical distribution object.
#' @param ... Additional arguments to pass (not used))
#' @return A function that takes a numeric vector `t` and returns the
#'         empirical cdf of `x` evaluated at `t`.
#' @importFrom stats ecdf
#' @export
cdf.empirical_dist <- function(x, ...) {
  if (dim(x) > 1) {
    stop("cdf not defined for multivariate empirical distribution")
  }
  ecdf(x$data)
}

#' Method for obtaining the support of `empirical_dist` object `x`.
#' @param x The empirical distribution object.
#' @return A `finite_set` object containing the support of `x`.
#' @export
sup.empirical_dist <- function(x) {
  finite_set$new(x$data)
}

#' Method for obtaining the number of observations used to construct a
#' `empirical_dist` object.
#' @param object The empirical distribution object.
#' @param ... Additional arguments to pass (not used).
#' @export
nobs.empirical_dist <- function(object, ...) {
  nrow(object$data)
}

#' Method for obtaining the observations used to construct a
#' `empirical_dist` object.
#' @param x The empirical distribution object.
#' @export
obs.empirical_dist <- function(x) {
  x$data
}


#' Method for obtaining the name of a `empirical_dist` object. Since the
#' empirical distribution is parameter-free, this function returns 0.
#' @param x The empirical distribution object.
#' @export
nparams.empirical_dist <- function(x) {
  0
}

#' `empirical_dist` objects have no parameters, so this function returns NULL.
#' @param x The empirical distribution object.
#' @export
params.empirical_dist <- function(x) {
  NULL
}

#' Method for printing a `dist` object
#' @param x The object to print
#' @param ... Additional arguments to pass
#' @export
print.empirical_dist <- function(x, ...) {
  print(summary(x, name = "Empirical Distribution", nobs = nobs(x), ...))
}



