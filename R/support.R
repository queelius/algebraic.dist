#' @title Support
#' @description support is a class that represents the support of a random
#' variable or distribution, i.e. the set of values that the random variable
#' can take on (non-zero probability or density).
#' 
#' It's a conceptual class, which a set of generic methods that can be
#' implemented for different types of support. To satisfy the concept of
#' a support, the following methods must be implemented:
#' 
#' 1. has: a function that returns a logical vector indicating
#' whether each value in a vector is contained in the support
#' 2. infimum: a function that returns the infimum of the support
#' 3. supremum: a function that returns the supremum of the support
#' 4. dim: a function that returns the dimension of the support
#' 
#' The following methods are also useful:
#' 
#' 1. print: a function that prints the support
#' 2. plot: a function that plots the support
#' 3. c: a function that combines two supports
#' 
#'We provide two implementations that satisfy the concept:
#' 
#' 1. `interval`: a support that is a finite union of intervals
#' 2. `infinite_set`: a support that is a finite set of values

#' Determine if a value is contained in the support.
#' @param object A support object.
#' @param x A vector of values.
#' @export
has <- function(object, x) {
    UseMethod("has", object)
}

#' Get the infimum of the support.
#' @param object A support object.
#' @export
infimum <- function(object) {
    UseMethod("infimum", object)
}

#' Get the supremum of the support.
#' @param object A support object.
#' @export
supremum <- function(object) {
    UseMethod("supremum", object)
}
