#' @title Interval
#' @description An interval is a support that is a finite union of intervals.
#' @field lower A numeric vector of lower bounds.
#' @field upper A numeric vector of upper bounds.
#' @field lower_closed A logical vector indicating whether the lower bound is
#' closed.
#' @field upper_closed A logical vector indicating whether the upper bound is
#' closed.
#' @importFrom R6 R6Class
#' @export
interval <- R6::R6Class(
    "interval",
    public = list(
        lower = NULL,
        upper = NULL,
        lower_closed = TRUE,
        upper_closed = TRUE,

        #' @description
        #' Initialize an interval.
        #' 
        #' @param lower A numeric vector of lower bounds.
        #' @param upper A numeric vector of upper bounds.
        #' @param lower_closed A logical vector indicating whether the lower
        #' bound is closed.
        #' @param upper_closed A logical vector indicating whether the upper
        #' bound is closed.
        initialize = function(lower = -Inf, upper = Inf,
                              lower_closed = FALSE, upper_closed = FALSE) {
            # replicate lower or upper if necessary so that they are the same
            # length
            if (length(lower) < length(upper)) {
                lower <- rep(lower, length(upper), length.out=length(upper))
            } else if (length(upper) < length(lower)) {
                upper <- rep(upper, length.out=length(lower))
            }
            # replicate lower_closed and upper_closed if necessary so that they
            # are the same length as `lower`
            if (length(lower_closed) != length(lower)) {
                lower_closed <- rep(lower_closed, length.out=length(lower))
            }
            if (length(upper_closed) != length(lower)) {
                upper_closed <- rep(upper_closed, length.out=length(lower))
            }

            self$lower <- lower
            self$upper <- upper
            self$lower_closed <- lower_closed
            self$upper_closed <- upper_closed
        },

        #' @description
        #' Determine if the interval is empty
        #' 
        #' @return A logical vector indicating whether the interval is empty.
        is_empty = function() {
            self$lower > self$upper || (
                self$lower == self$upper &
                !(self$lower_closed & self$upper_closed))
        },

        #' @description
        #' Determine if a value is contained in the interval.
        #' 
        #' @param x A numeric vector of values.
        #' @return A logical vector indicating whether each value is contained
        contains = function(x) {
            lower <- ifelse(self$lower_closed, x >= self$lower, x > self$lower)
            upper <- ifelse(self$upper_closed, x <= self$upper, x < self$upper)
            lower & upper
        },

        #' @description
        #' Get the infimum of the interval.
        #' 
        #' @return A numeric vector of infimums.
        infimum = function() {
            self$lower
        },

        #' @description
        #' Get the supremum of the interval.
        #' 
        #' @return A numeric vector of supremums.
        supremum = function() {
            self$upper
        },

        #' @description
        #' Get the dimension of the interval.
        #' 
        #' @return The dimension of the interval.
        dim = function() {
            length(self$lower)
        }
    )
)

#' Determine if a value is contained in the interval.
#' @param object An interval object.
#' @param x A vector of values.
#' @export
contains.interval <- function(object, x) {
    object$contains(x)
}

#' Return the (vector of) infimum of the interval.
#' @param object An interval object.
#' @export
infimum.interval <- function(object) {
    object$infimum()
}

#' Return the (vector of) supremum of the interval.
#' @param object An interval object.
#' @export
supremum.interval <- function(object) {
    object$supremum()
}

#' Return the dimension of the interval.
#' @param object An interval object.
#' @export
dim.interval <- function(object) {
    object$dim()
}

#' Print the interval.
#' @param x An interval object.
#' @param ... Additional arguments.
#' @export
print.interval <- function(x, ...) {
    for (i in 1:length(x$lower)) {
        if (x$lower_closed[i]) cat("[", sep = "")
        else cat("(", sep = "")
        cat(x$lower[i], ", ", x$upper[i], sep = "")
        if (x$upper_closed[i]) cat("]", sep = "")
        else cat(")", sep = "")
        cat("\n")
    }
}
