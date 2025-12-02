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
#' Attempts to reduce expression distributions to closed-form distributions
#' when mathematical identities apply. For example:
#' - normal + normal = normal (sum of independent normals)
#' - normal - normal = normal (difference of independent normals)
#'
#' @param x The `edist` object to simplify
#' @param ... Additional arguments to pass (not used)
#' @return The simplified distribution, or unchanged `edist` if no rule applies
#' @export
simplify.edist <- function(x, ...) {
  expr <- x$e
  vars <- x$vars

  # Check for binary operations with two operands

  if (length(vars) == 2) {
    d1 <- vars[[1]]
    d2 <- vars[[2]]

    # Detect operation from expression
    op <- if (is.call(expr)) as.character(expr[[1]]) else NULL

    # Rule: normal + normal -> normal
    # X ~ N(μ₁, σ₁²), Y ~ N(μ₂, σ₂²) => X + Y ~ N(μ₁ + μ₂, σ₁² + σ₂²)
    if (identical(op, "+") && is_normal(d1) && is_normal(d2)) {
      return(normal(
        mu = mean(d1) + mean(d2),
        var = vcov(d1) + vcov(d2)
      ))
    }

    # Rule: normal - normal -> normal
    # X ~ N(μ₁, σ₁²), Y ~ N(μ₂, σ₂²) => X - Y ~ N(μ₁ - μ₂, σ₁² + σ₂²)
    if (identical(op, "-") && is_normal(d1) && is_normal(d2)) {
      return(normal(
        mu = mean(d1) - mean(d2),
        var = vcov(d1) + vcov(d2)
      ))
    }

    # Additional rules can be added here:
    # - exponential + exponential (same rate) -> Erlang/Gamma
    # - sum of n iid exponentials -> Gamma(n, rate)
    # - etc.
  }

  # No simplification rule matched, return unchanged
  x
}



#' We have an edist object, which is a subclass of dist, and now we're placing
#' it inside of a limit expression, where the limit is understood to be with
#' respect to sample size. We need to define a method for this.
#'
#' @param x The edist object to take the limit of
#' @return The limit of the edist object
#' @keywords internal
limit.edist <- function(x) {
  # we just wrap the edist object in a limit object, and the simplify method
  # can be used to simplify the limit object if necessary

  expression(limit(x$e))
}