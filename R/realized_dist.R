#' Construct a realized distribution object.
#'
#' A \code{realized_dist} is an \code{\link{empirical_dist}} that preserves
#' provenance: it remembers the source distribution that generated its
#' samples.  All \code{empirical_dist} methods work unchanged via S3
#' inheritance.
#'
#' @param data A matrix of sampled data (rows = observations).
#' @param source The distribution object that generated \code{data}.
#' @param n The number of samples drawn.
#' @return A \code{realized_dist} object.
#' @keywords internal
realized_dist <- function(data, source, n) {
  ed <- empirical_dist(data)
  ed$source <- source
  ed$n_realized <- n
  class(ed) <- c("realized_dist", class(ed))
  ed
}

#' Test whether an object is a \code{realized_dist}.
#'
#' @param x The object to test.
#' @return \code{TRUE} if \code{x} inherits from \code{"realized_dist"},
#'   \code{FALSE} otherwise.
#' @export
is_realized_dist <- function(x) {
  inherits(x, "realized_dist")
}

#' Format a \code{realized_dist} object as a character string.
#'
#' Shows the number of samples and a summary of the source distribution.
#'
#' @param x A \code{realized_dist} object.
#' @param ... Additional arguments (not used).
#' @return A character string.
#' @export
format.realized_dist <- function(x, ...) {
  sprintf("Realized distribution (%d samples from: %s)",
          x$n_realized, format(x$source))
}

#' Print a \code{realized_dist} object.
#'
#' @param x A \code{realized_dist} object.
#' @param ... Additional arguments (not used).
#' @return \code{x}, invisibly.
#' @export
print.realized_dist <- function(x, ...) {
  cat(format(x), "\n")
  cat("  source:\n")
  print(x$source)
  invisible(x)
}
