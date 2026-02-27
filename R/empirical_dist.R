#' Construct empirical distribution object.
#' @param data data to construct empirical distribution from. if matrix or
#'             data frame, each row is a joint observation, if a vector, each
#'             element is an observation. whatever data is, it must be
#'             convertible to a tibble.
#' @return An \code{empirical_dist} object.
#' @examples
#' # Univariate empirical distribution from a vector
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' mean(ed)
#'
#' # Multivariate empirical distribution from a matrix
#' mat <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' ed_mv <- empirical_dist(mat)
#' dim(ed_mv)
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
#' @return Logical; \code{TRUE} if \code{x} is an \code{empirical_dist}.
#' @examples
#' ed <- empirical_dist(c(1, 2, 3))
#' is_empirical_dist(ed)    # TRUE
#' is_empirical_dist("abc") # FALSE
#' @export
is_empirical_dist <- function(x) {
  inherits(x, "empirical_dist")
}

#' Method for obtaining the dimension of a `empirical_dist` object.
#' @param x The object to obtain the dimension of.
#' @return Integer; the number of dimensions.
#' @examples
#' ed1 <- empirical_dist(c(1, 2, 3))
#' dim(ed1) # 1
#'
#' ed2 <- empirical_dist(matrix(1:6, ncol = 2))
#' dim(ed2) # 2
#' @export
dim.empirical_dist <- function(x) {
  ncol(x$data)
}

#' Method for obtaining the pdf of a `empirical_dist` object.
#'
#' @param x The object to obtain the pdf of.
#' @param ... Additional arguments to pass into the pdf function.
#' @note sort tibble lexicographically and do a binary search to find upper
#'       and lower bound in `log(nobs(x))` time.
#' @return A function computing the empirical PMF at given points.
#' @examples
#' ed <- empirical_dist(c(1, 2, 2, 3, 3, 3))
#' f <- density(ed)
#' f(2)          # 2/6
#' f(3, log = TRUE) # log(3/6)
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
#' @return A function that takes \code{n} and returns \code{n} resampled
#'   observations.
#' @examples
#' ed <- empirical_dist(c(10, 20, 30))
#' s <- sampler(ed)
#' set.seed(42)
#' s(5)
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
#' @param ... Additional arguments to pass into function `g`.
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
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' expectation(ed)                     # E[X] = 3
#' expectation(ed, function(x) x^2)    # E[X^2] = 11
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
#' @return Numeric vector of column means.
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' mean(ed) # 3
#' @export
mean.empirical_dist <- function(x, ...) {
  colMeans(x$data)
}

#' Method for obtaining the variance of `empirical_dist` object `x`.
#'
#' @param object The empirical distribution object.
#' @param ... Additional arguments to pass (not used).
#' @return The sample variance-covariance matrix.
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' vcov(ed) # sample variance
#'
#' ed_mv <- empirical_dist(matrix(rnorm(20), ncol = 2))
#' vcov(ed_mv) # 2x2 covariance matrix
#' @importFrom stats cov
#' @export
vcov.empirical_dist <- function(object, ...) {
  cov(object$data)
}

#' Method for obtaining the marginal distribution of `empirical_dist` object
#' `x`.
#' @param x The empirical distribution object.
#' @param indices The indices of the marginal distribution to obtain.
#' @return An \code{empirical_dist} over the selected columns.
#' @examples
#' mat <- matrix(1:12, ncol = 3)
#' ed <- empirical_dist(mat)
#' ed_marginal <- marginal(ed, c(1, 3))
#' dim(ed_marginal) # 2
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
#' @return An \code{empirical_dist} containing only rows satisfying \code{P}.
#' @examples
#' \donttest{
#' mat <- matrix(c(1, 5, 2, 3, 4, 1, 6, 2), ncol = 2)
#' ed <- empirical_dist(mat)
#' # Condition on first column being greater than second
#' ed_cond <- conditional(ed, function(d) d[1] > d[2])
#' nobs(ed_cond)
#' }
#' @export
conditional.empirical_dist <- function(x, P, ...) {
  mask <- apply(X = x$data, MARGIN = 1, FUN = P, ...)
  if (!any(mask))
    stop("conditioning resulted in zero observations; predicate matched no rows")
  empirical_dist(x$data[mask, , drop = FALSE])
}


#' Method for obtaining the empirical distribution of a function of the
#' observations of `empirical_dist` object `x`.
#' 
#' @param x The empirical distribution object.
#' @param g The function to apply to each observation.
#' @param ... Additional arguments to pass into function `g`.
#' @return An \code{empirical_dist} of the transformed observations.
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4))
#' ed2 <- rmap(ed, function(x) x^2)
#' mean(ed2) # mean of 1, 4, 9, 16
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
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' Fx <- cdf(ed)
#' Fx(3) # 0.6
#' Fx(c(1, 5)) # 0.2, 1.0
#' @importFrom stats ecdf
#' @export
cdf.empirical_dist <- function(x, ...) {
  if (dim(x) > 1) {
    stop("cdf not defined for multivariate empirical distribution")
  }
  ecdf(x$data)
}

#' Method for obtaining the inverse CDF (quantile function) of a univariate
#' `empirical_dist` object.
#'
#' Uses the empirical quantile function from the observed data.
#'
#' @param x The empirical distribution object.
#' @param ... Additional arguments (not used).
#' @return A function that accepts a vector of probabilities `p` and returns
#'   the corresponding quantiles.
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' qf <- inv_cdf(ed)
#' qf(0.5)       # median
#' qf(c(0.25, 0.75)) # quartiles
#' @importFrom stats quantile
#' @export
inv_cdf.empirical_dist <- function(x, ...) {
  if (dim(x) > 1) {
    stop("inv_cdf not defined for multivariate empirical distribution")
  }
  xdata <- as.numeric(x$data)
  function(p, ...) {
    as.numeric(quantile(xdata, probs = p, type = 1, ...))
  }
}

#' Method for obtaining the support of `empirical_dist` object `x`.
#' @param x The empirical distribution object.
#' @return A `finite_set` object containing the support of `x`.
#' @examples
#' ed <- empirical_dist(c(1, 2, 2, 3))
#' s <- sup(ed)
#' s$has(2) # TRUE
#' s$has(4) # FALSE
#' @export
sup.empirical_dist <- function(x) {
  finite_set$new(x$data)
}

#' Method for obtaining the number of observations used to construct a
#' `empirical_dist` object.
#' @param object The empirical distribution object.
#' @param ... Additional arguments to pass (not used).
#' @return Integer; number of observations.
#' @examples
#' ed <- empirical_dist(c(10, 20, 30, 40))
#' nobs(ed) # 4
#' @export
nobs.empirical_dist <- function(object, ...) {
  nrow(object$data)
}

#' Method for obtaining the observations used to construct a
#' `empirical_dist` object.
#' @param x The empirical distribution object.
#' @return A matrix of observations (rows = observations, columns = dimensions).
#' @examples
#' ed <- empirical_dist(c(5, 10, 15))
#' obs(ed)
#' @export
obs.empirical_dist <- function(x) {
  x$data
}


#' Method for obtaining the name of a `empirical_dist` object. Since the
#' empirical distribution is parameter-free, this function returns 0.
#' @param x The empirical distribution object.
#' @return 0 (empirical distributions are non-parametric).
#' @examples
#' ed <- empirical_dist(c(1, 2, 3))
#' nparams(ed) # 0
#' @export
nparams.empirical_dist <- function(x) {
  0
}

#' `empirical_dist` objects have no parameters, so this function returns NULL.
#' @param x The empirical distribution object.
#' @return \code{NULL} (empirical distributions have no parameters).
#' @examples
#' ed <- empirical_dist(c(1, 2, 3))
#' params(ed) # NULL
#' @export
params.empirical_dist <- function(x) {
  NULL
}

#' Format method for `empirical_dist` objects.
#' @param x The object to format
#' @param ... Additional arguments (not used)
#' @return A character string
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' format(ed)
#' @export
format.empirical_dist <- function(x, ...) {
  sprintf("Empirical distribution (%d observations, %d dimensions)",
          nobs(x), ncol(x$data))
}

#' Print method for `empirical_dist` objects.
#' @param x The object to print
#' @param ... Additional arguments to pass
#' @return \code{x}, invisibly.
#' @examples
#' ed <- empirical_dist(c(1, 2, 3, 4, 5))
#' print(ed)
#' @export
print.empirical_dist <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}



