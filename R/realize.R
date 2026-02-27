#' Materialize any distribution to empirical_dist by sampling.
#'
#' \code{realize} draws \code{n} samples from a distribution and wraps
#' them in an \code{\link{empirical_dist}}.
#' This is the universal fallback that lets any \code{dist} object be
#' converted to a discrete approximation on which methods like
#' \code{\link{cdf}}, \code{\link{density}}, and \code{\link{conditional}}
#' are always available.
#'
#' For non-empirical distributions, the result is a
#' \code{\link{realized_dist}} that preserves the source distribution
#' as provenance metadata.  This enables re-sampling via
#' \code{realize(x$source, n = ...)} and informative printing.
#'
#' The \code{empirical_dist} method is a no-op: the distribution is
#' already materialized.
#'
#' The \code{realized_dist} method re-samples from the original source
#' distribution, allowing cheap regeneration with a different sample size.
#'
#' @param x A distribution object (inheriting from \code{dist}).
#' @param n Number of samples (default: 10000).
#' @param ... Additional arguments passed to methods.
#' @return An \code{\link{empirical_dist}} (or \code{\link{realized_dist}})
#'   object.
#' @export
realize <- function(x, n = 10000, ...) UseMethod("realize")

#' @rdname realize
#' @export
realize.dist <- function(x, n = 10000, ...) {
  n <- validate_n(n)
  samples <- sampler(x)(n)
  realized_dist(samples, source = x, n = n)
}

#' @rdname realize
#' @export
realize.empirical_dist <- function(x, ...) {
  x
}

#' @rdname realize
#' @export
realize.realized_dist <- function(x, n = 10000, ...) {
  n <- validate_n(n)
  realize(x$source, n = n)
}

validate_n <- function(n) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("'n' must be a positive integer, got: ", deparse(n))
  as.integer(n)
}


# ---- Internal memoized fallback entry point --------------------------------

#' Memoized MC fallback materialization.
#'
#' Single internal entry point for all Monte Carlo fallback paths.
#' If \code{x} is already an \code{empirical_dist}, returns it unchanged.
#' If \code{x} has a \code{.cache} environment (e.g. \code{edist} objects),
#' caches the realization so that multiple method calls (e.g. \code{cdf} +
#' \code{density}) share the same samples.  Sample-size-aware: if the cached
#' realization has fewer than \code{n} samples, re-realizes.
#'
#' @param x A distribution object.
#' @param n Number of samples (default: 10000).
#' @return An \code{empirical_dist} (or \code{realized_dist}).
#' @keywords internal
ensure_realized <- function(x, n = 10000L) {
  if (is_empirical_dist(x)) return(x)

  cache <- x$.cache
  if (is.environment(cache)) {
    cached <- cache$.realized
    if (!is.null(cached) && cached$n_realized >= n) {
      return(cached)
    }
    result <- realize(x, n = n)
    cache$.realized <- result
    return(result)
  }

  realize(x, n = n)
}
