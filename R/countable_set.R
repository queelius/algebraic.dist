#' @title Countable Set
#' @description A countably infinite support set, such as the non-negative
#' integers. It satisfies the concept of a support (see \code{\link{has}},
#' \code{\link{infimum}}, \code{\link{supremum}}, \code{\link[base]{dim}}).
#'
#' @field lower_bound Integer lower bound of the set.
#' @importFrom R6 R6Class
#' @export
countable_set <- R6::R6Class("countable_set",
  public = list(
    lower_bound = NULL,

    #' @description
    #' Initialize a countable set.
    #'
    #' @param lower Integer lower bound (default 0).
    initialize = function(lower = 0L) {
      self$lower_bound <- as.integer(lower)
    }
  )
)

#' Check membership in a countable set.
#'
#' Returns \code{TRUE} if all values are integers (within floating-point
#' tolerance) that are at least as large as the lower bound.
#'
#' @param object A \code{countable_set} object.
#' @param x Value(s) to check.
#' @return Logical; \code{TRUE} if all values are valid members of the set.
#' @examples
#' cs <- countable_set$new(0L)
#' has(cs, c(0, 3, 5))   # TRUE
#' has(cs, c(-1, 2))     # FALSE (negative integer)
#' has(cs, 1.5)          # FALSE (not integer)
#' @export
has.countable_set <- function(object, x) {
  all(is.numeric(x) & x >= object$lower_bound & x == floor(x))
}

#' Get the infimum of a countable set.
#'
#' @param object A \code{countable_set} object.
#' @return The lower bound (integer).
#' @examples
#' cs <- countable_set$new(0L)
#' infimum(cs)  # 0
#' @export
infimum.countable_set <- function(object) {
  object$lower_bound
}

#' Get the supremum of a countable set.
#'
#' @param object A \code{countable_set} object.
#' @return \code{Inf} (the set is unbounded above).
#' @examples
#' cs <- countable_set$new(0L)
#' supremum(cs)  # Inf
#' @export
supremum.countable_set <- function(object) {
  Inf
}

#' Get the dimension of a countable set.
#'
#' @param x A \code{countable_set} object.
#' @return \code{1} (always univariate).
#' @examples
#' cs <- countable_set$new(0L)
#' dim(cs)  # 1
#' @export
dim.countable_set <- function(x) {
  1
}
