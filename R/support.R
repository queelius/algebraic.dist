#' @export
contains <- function(support, x) {
    UseMethod("contains", support)
}

#' @export
infimum <- function(support) {
    UseMethod("infimum", support)
}

#' @export
supremum <- function(support) {
    UseMethod("supremum", support)
}

interval <- R6::R6Class(
    classname = c("interval", "set"),
    public = list(
        lower = NULL,
        upper = NULL,
        lower_closed = TRUE,
        upper_closed = TRUE,
        initialize = function(lower = -Inf, upper = Inf,
                              lower_closed = TRUE, upper_closed = TRUE) {
            # replicate lower or upper if necessary so that they are the same length
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
        contains = function(x) {
            lower <- ifelse(self$lower_closed, x >= self$lower, x > self$lower)
            upper <- ifelse(self$upper_closed, x <= self$upper, x < self$upper)
            lower & upper
        },
        infimum = function() {
            self$lower
        },
        supremum = function() {
            self$upper
        },
        dim = function() {
            length(self$lower)
        }
    )
)

#' @export
contains.interval <- function(support, x) {
    support$contains(x)
}

#' @export
infimum.interval <- function(support) {
    support$infimum()
}

#' @export
supremum.interval <- function(support) {
    support$supremum()
}

#' @export
dim.interval <- function(support) {
    support$dim()
}

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

finite_set <- R6::R6Class(
    classname = c("finite_set", "set"),
    public = list(
        values = NULL,
        initialize = function(values) {
            self$values <- unique(values)
        },
        contains = function(x) {
            if (is.matrix(self$values)) {
                # check if x is the same as any row in values matrix
                return(any(apply(self$values, 1, function(row) {
                    all(row == x)
                })))
            } else {
                return(x %in% self$values)
            }
        },
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
        dim = function() {
            if (is.matrix(self$values)) {
                ncol(self$values)
            } else {
                1
            }
        }
    )
)

#' @export
contains.finite_set <- function(support, x) {
    support$contains(x)
}

#' @export
infimum.finite_set <- function(support) {
    support$infimum()
}

#' @export
supremum.finite_set <- function(support) {
    support$supremum()
}

#' @export
dim.finite_set <- function(support) {
    support$dim()
}

