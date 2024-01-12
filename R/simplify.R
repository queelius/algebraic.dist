#' Generic method for simplifying distributions.
#' 
#' @param x The distribution to simplify
#' @return The simplified distribution
#' @export
simplify <- function(x) UseMethod("simplify")


#' Default Method for simplifming a `dist` object.
#' @param x The `dist` object to simplify
#' @return The `dist` object
#' @export
simplify.dist <- function(x) {
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
#' @return The simplified object
#' @export
#' @examples
#' e <- edist(expression(x + y),
#'            list(x = exponential(rate = 2),
#'                 y = exponential(rate = 1)))
#' print(e)
#' #> Distribution x + y
#' #>   x ~ exponential(rate = 2)
#' #>   y ~ exponential(rate = 1)
#' e.simp <- simplify(e)
#' print(e.simp)
#' #> Exponential distribution with rate 3
simplify.edist <- function(x) {

    # First, we need to walk the expression tree from the bottom up
    # and try to simplify each node. We do this by recursively calling
    # simplify on each node of the expression tree.
    x <- lapply(x, simplify)
    
    # Next, we need to try all valid simplifications and choose the
    # one that is the most simplified. We do this by creating a list
    # of all valid simplifications and then using the `simplify` function
    # to choose the most simplified one.
    simplifications <- list(
        simplify.min.exponential,
        simplify.sum.exponential,
        simplify.constant
    )
    x <- simplify(x, simplifications)
    
    # Finally, we need to update the class of the distribution.
    # We do this by calling the `update_class` function.
    update_class(x)

}
