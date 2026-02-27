#' @title Finite set
#' @description A finite set. It also satisfies the concept of a support.
#' 
#' @field values A vector of values.
#' @importFrom R6 R6Class
#' @export
finite_set <- R6::R6Class(
    "finite_set",
    public = list(
        values = NULL,

        #' @description
        #' Initialize a finite set.
        #' 
        #' @param values A vector of values.
        initialize = function(values) {
            self$values <- unique(values)
        },

        #' @description
        #' Determine if a value is contained in the finite set.
        #' 
        #' @param x A vector of values.
        has = function(x) {
            if (is.matrix(self$values)) {
                # check if x is the same as any row in values matrix
                return(any(apply(self$values, 1, function(row) {
                    all(row == x)
                })))
            } else {
                return(x %in% self$values)
            }
        },

        #' @description
        #' Get the infimum of the finite set.
        #' 
        #' @return A numeric vector of infimums.
        infimum = function() {
             # if values is a matrix, find the min of each column
            if (is.matrix(self$values)) {
                sapply(seq_len(ncol(self$values)), function(i) {
                    min(self$values[, i])
                })
            } else {
                min(self$values)
            }
        },

        #' @description
        #' Get the supremum of the finite set.
        #' 
        #' @return A numeric vector of supremums.
        supremum = function() {
            # if values is a matrix, find the max of each column
            if (is.matrix(self$values)) {
                sapply(seq_len(ncol(self$values)), function(i) {
                    max(self$values[, i])
                })
            } else {
                max(self$values)
            }
        },

        #' @description
        #' Get the dimension of the finite set.
        #' 
        #' @return The dimension of the finite set.
        dim = function() {
            if (is.matrix(self$values)) {
                ncol(self$values)
            } else {
                1
            }
        }
    )
)

#' Determine if a value is contained in the finite set.
#' @param object A finite set.
#' @param x A vector of values.
#' @return Logical indicating membership.
#' @examples
#' fs <- finite_set$new(c(1, 3, 5, 7))
#' has(fs, 3) # TRUE
#' has(fs, 4) # FALSE
#' @export
has.finite_set <- function(object, x) {
    object$has(x)
}

#' Return the infimum of the finite set.
#' @param object A finite set.
#' @return Numeric; the minimum value(s).
#' @examples
#' fs <- finite_set$new(c(1, 3, 5, 7))
#' infimum(fs) # 1
#' @export
infimum.finite_set <- function(object) {
    object$infimum()
}

#' Return the supremum of the finite set.
#' @param object A finite set.
#' @return Numeric; the maximum value(s).
#' @examples
#' fs <- finite_set$new(c(1, 3, 5, 7))
#' supremum(fs) # 7
#' @export
supremum.finite_set <- function(object) {
    object$supremum()
}

#' Return the dimension of the finite set.
#' @param x A finite set.
#' @return Integer; the dimension of the set.
#' @examples
#' fs <- finite_set$new(c(1, 3, 5, 7))
#' dim(fs) # 1
#' @export
dim.finite_set <- function(x) {
    x$dim()
}

