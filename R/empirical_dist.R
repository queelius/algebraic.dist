#' Construct empirical distribution object.
#' 
#' @param data data to construct empirical distribution from. if matrix or
#'             data frame, each row is a joint observation, if a vector, each
#'             element is an observation. whatever data is, it must be
#'             convertible to a tibble.
#' @export
empirical_dist <- function(data) {

  cls_names <- if (is.matrix(data)) {
    c("empirical_dist", "multivariate_dist", "dist")
  } else {
    c("empirical_dist", "univariate_dist", "dist")
  }
  structure(list(data = as.matrix(data)), class = cls_names)
}

#' Function to determine whether an object `x` is an `empirical_dist` object.
#' @export
is_empirical_dist <- function(x) {
  inherits(x, "empirical_dist")
}

#' Method for obtaining the number of observations in a `empirical_dist`
#' object.
#' @param x The object to obtain the number of observations of.
#' @export
nobs.empirical_dist <- function(x) {
  nrow(x$data)
}


#' Method for obtaining the dimension of a `empirical_dist` object.
#' @param x The object to obtain the dimension of.
#' @export
dim.empirical_dist <- function(x) {
  ncol(x$data)
}

#' Method for obtaining the data of a `empirical_dist` object.
#' @param x The object to obtain the data of.
#' @export
obs.empirical_dist <- function(x) {
  x$data
}

#' Method for obtaining the pdf of a `empirical_dist` object.
#'
#' @param x The object to obtain the pdf of.
#' @note sort tibble lexographically and do a binary search to find upper
#'       and lower bound in `log(nobs(x))` time.
#' @export
pdf.empirical_dist <- function(x) {
  n <- nrow(x$data)
  p <- ncol(x$data)
  function(t, log = FALSE) {
    stopifnot(is.numeric(t), length(t) == p)

    # return number of elements in x equal to t
    count <- sum(x$data == t)
    if (log) {
      return(log(count) - log(n))
    } else {
      return(count / n)
    }
  }
}

#' Method for obtaining the sampler for a `empirical_dist` object.
#'
#' @param x The object to obtain the sampler of.
#' @export
sampler.empirical_dist <- function(x) {
  N <- nrow(x$data)
  function (n = 1, ...) {
    x$data[sample(1:N, size = n, replace = TRUE, ...), ]
  }
}

#' Method for obtaining the expectation of `empirical_dist` object `x`
#' under function `g`.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass itno function `g`.
#' @export
expectation.empirical_dist <- function(x, g, ...) {
  # apply g to each row of x
  # and take their column means
  colMeans(apply(x$data, 1, g, ...))
}

#' Method for obtaining the mean of `empirical_dist` object `x`.
#'
#' @param x The distribution object.
#' @export
mean.empirical_dist <- function(x) {
  colMeans(x$data)
}

#' Method for obtaining the variance of `empirical_dist` object `x`.
#'
#' @param x The empirical distribution object.
#' @export
vcov.empirical_dist <- function(x) {
  cov(x$data)
}

#' Method for obtaining the marginal distribution of `empirical_dist` object
#' `x`.
#' @param x The empirical distribution object.
#' @param i The indices of the marginal distribution to obtain.
#' @export
marginal.empirical_dist <- function(x, indices) {
  empirical_dist(x$data[, indices])
}

#' Method for obtaining the condition distribution, x | P(x), of
#' `empirical_dist` object `x`.
#' 
#' In other words, we condition the data on the predicate function `P`.
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
#' @param P the predicate function we condition the data on.
conditional.empirical_dist <- function(x, P) {
  x$data <- x$data[apply(x$data, 1, P), ]
  x
}


#' Method for obtaining the empirical distribution of a function of the
#' observations of `empirical_dist` object `x`.
#' 
#' @param x The empirical distribution object.
#' @param g The function to apply to each observation.
#' @export
rmap.empirical_dist <- function(x, g) {
  x$data <- apply(x$data, 1, g)
  x
}

#' Method for obtaining the cdf of `empirical_dist` object `x`.
#'
#' If `x` is a multivariate empirical distribution, this function will
#' throw an error. It's only defined for univariate empirical distributions.
#' @param x The empirical distribution object.
#' @return A function that takes a numeric vector `t` and returns the
#'         empirical cdf of `x` evaluated at `t`.
#' @export
cdf.empirical_dist <- function(x) {
  if (dim(x) > 1) {
    stop("cdf not defined for multivariate empirical distribution")
  }
  ecdf(x$data)
}

#' @export
sup.empirical_dist <- function(x) {
  finite_set$new(x$data)
}