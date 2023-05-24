#' @export
contains <- function(support, x) {
    UseMethod("contains", support)
}

#' @export
infinum <- function(support) {
    UseMethod("infinum", support)
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
            self$lower <- lower
            self$upper <- upper
            self$lower_closed <- lower_closed
            self$upper_closed <- upper_closed
        },
        contains = function(x) {
            if (self$lower_closed && self$upper_closed) {
                return(self$lower <= x & x <= self$upper)
            } else if (self$lower_closed) {
                return(self$lower <= x & x < self$upper)
            } else if (self$upper_closed) {
                return(self$lower < x & x <= self$upper)
            } else {
                return(self$lower < x & x < self$upper)
            }
        },
        infinum = function() {
            self$lower
        },
        supremum = function() {
            self$upper
        },
        dim = function() {
            1
        }
    )
)

#' @export
contains.interval <- function(support, x) {
    support$contains(x)
}

#' @export
infinum.interval <- function(support) {
    support$infinum()
}

#' @export
supremum.interval <- function(support) {
    support$supremum()
}

#' @export
dim.interval <- function(support) {
    support$dim()
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
        infinum = function() {
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
infinum.finite_set <- function(support) {
    support$infinum()
}

#' @export
supremum.finite_set <- function(support) {
    support$supremum()
}

#' @export
dim.finite_set <- function(support) {
    support$dim()
}

#' @export
box_support <- function(dim = 1, lower = -Inf, upper = Inf) {
    structure(list(lower = lower, upper = upper, dim = dim),
             class = c("box_support", "set"))
}

#' @export
contains.box_support <- function(support, x) {
    all(x >= support$lower & x <= support$upper)
}

#' @export
infinum.box_support <- function(support) {
    rep(support$lower, support$dim)
}

#' @export
supremum.box_support <- function(support) {
    rep(support$upper, support$dim)
}

#' @export
dim.box_support <- function(support) {
    support$dim
}

