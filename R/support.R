#' @title Support
#' @description support is a class that represents the support of a random
#' element or distribution, i.e. the set of values that it realize.
#' 
#' It's a conceptual class. To satisfy the concept of
#' a support, the following methods must be implemented:
#' 
#' 1. has: a function that returns a logical vector indicating
#' whether each value in a vector is contained in the support
#' 2. infimum: a function that returns the infimum of the support
#' 3. supremum: a function that returns the supremum of the support
#' 4. dim: a function that returns the dimension of the support
#' 
#'We provide two implementations that satisfy the concept:
#' 
#' - `interval`: a support that is an infiite set of contiguous numeric values
#' - `finite_set`: a support that is a finite set of values

#' Determine if a value is contained in the support.
#' @param object A support object.
#' @param x A vector of values.
#' @return Logical vector indicating membership.
#' @examples
#' I <- interval$new(0, 1, lower_closed = TRUE, upper_closed = TRUE)
#' has(I, 0.5)  # TRUE
#' has(I, 2)    # FALSE
#'
#' S <- finite_set$new(c(1, 2, 3))
#' has(S, 2)    # TRUE
#' has(S, 4)    # FALSE
#' @export
has <- function(object, x) {
    UseMethod("has", object)
}

#' Get the infimum of the support.
#' @param object A support object.
#' @return The infimum (greatest lower bound) of the support.
#' @examples
#' I <- interval$new(0, 10)
#' infimum(I)  # 0
#'
#' S <- finite_set$new(c(3, 7, 11))
#' infimum(S)  # 3
#' @export
infimum <- function(object) {
    UseMethod("infimum", object)
}

#' Get the supremum of the support.
#' @param object A support object.
#' @return The supremum (least upper bound) of the support.
#' @examples
#' I <- interval$new(0, 10)
#' supremum(I)  # 10
#'
#' S <- finite_set$new(c(3, 7, 11))
#' supremum(S)  # 11
#' @export
supremum <- function(object) {
    UseMethod("supremum", object)
}
