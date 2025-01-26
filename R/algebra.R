#' Generic method for simplifying distributions.
#' 
#' @param x The distribution to simplify
#' @param ... Additional arguments to pass
#' @return The simplified distribution
#' @export
simplify <- function(x, ...) UseMethod("simplify")

#' Default Method for simplifming a `dist` object. Just returns the object.
#' @param x The `dist` object to simplify
#' @param ... Additional arguments to pass (not used)
#' @return The `dist` object
#' @export
simplify.dist <- function(x, ...) {
  x
}

#' Method for simplifying an `edist` object.
#' 
#' This is a complicated function that walks the expression tree
#' and tries to simplify it. Since sometimes a simplification
#' made at some level of the expression tree can lead to a
#' simplification at a higher level of the expression tree,
#' we need to walk the expression tree from the bottom up.
#' 
#' Also, since some simplifications can lead to a change in
#' the class of the distribution, we need to be careful to
#' update the class of the distribution as we simplify it.
#' 
#' Finally, the simplifications we initially choose to do can
#' prevent us from doing other simplifications that may ultimately
#' be more beneficial. So, we need to try all valid simplifications,
#' creating new `edist` objects for each simplification, and then
#' choose the one that is the most simplified.
#' 
#' @param x The `edist` object to simplify
#' @param ... Additional arguments to pass (not used)
#' @return The simplified object
#' @export
simplify.edist <- function(x, ...) {
   # we need to walk the expression tree from the bottom up
   # and try all valid simplifications
   # we also need to update the class of the distribution
   # as we simplify it
   # finally, we need to choose the most simplified object
   # from all the simplifications we tried
   # for now, we just return the object
   x
}



#' We have an edist object, which is a subclass of dist, and now we're placing
#' it inside of a limit expression, where the limit is understood to be with
#' respect to sample size. We need to define a method for this.
#' 
#' @param x The edist object to take the limit of
#' @param ... Additional arguments to pass (not used)
#' @return The limit of the edist object
limit.edist <- function(x) {
  # we just wrap the edist object in a limit object, and the simplify method
  # can be used to simplify the limit object if necessary

  expression(limit(x$e))
}