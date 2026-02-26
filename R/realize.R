#' Materialize any distribution to empirical_dist by sampling.
#'
#' \code{realize} draws \code{n} samples from a distribution and wraps
#' them in an \code{\link{empirical_dist}}.
#' This is the universal fallback that lets any \code{dist} object be
#' converted to a discrete approximation on which methods like
#' \code{\link{cdf}}, \code{\link{density}}, and \code{\link{conditional}}
#' are always available.
#'
#' The \code{empirical_dist} method is a no-op: the distribution is
#' already materialized.
#'
#' @param x A distribution object (inheriting from \code{dist}).
#' @param n Number of samples (default: 10000).
#' @param ... Additional arguments passed to methods.
#' @return An \code{\link{empirical_dist}} object.
#' @export
realize <- function(x, n = 10000, ...) UseMethod("realize")

#' @rdname realize
#' @export
realize.dist <- function(x, n = 10000, ...) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("'n' must be a positive integer, got: ", deparse(n))
  n <- as.integer(n)
  samples <- sampler(x)(n)
  empirical_dist(samples)
}

#' @rdname realize
#' @export
realize.empirical_dist <- function(x, ...) {
  x
}
