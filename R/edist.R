#' Takes an expression `e` and a list `vars` and returns a
#' lazy `edist` (expression distribution object), that is a subclass
#' of `dist` that can be used in place of a `dist` object.
#' 
#' @param e the expression to evaluate against the arguments.
#' @param vars the list of distributions (with variable names)
#'             to evaluate the expression `e` against.
#' @return An `edist` object.
#' @export
edist <- function(e, vars) {
  # retrieve the class of each of the distributions
  # to include that in the class of the edist object
  # we are creating
  classes <- sapply(vars, class)[1,]
  classes <- paste0(classes, collapse = "_")
  expr_str <- paste0(e, "_", classes)
  expr_str <- gsub(" ", "_", expr_str)
  structure(list(e = e, vars = vars, .cache = new.env(parent = emptyenv())),
            class = c(expr_str, "edist", "dist"))
}

#' Function to determine whether an object `x` is an `edist` object.
#' @param x The object to test
#' @export
is_edist <- function(x) {
  inherits(x, "edist")
}

#' Method for obtaining the parameters of an `edist` object.
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @export
params.edist <- function(x) {
  do.call(c, lapply(x$vars, params))
}

#' Method for obtaining the variance-covariance matrix (or scalar)
#' 
#' @param object The `edist` object to retrieve the variance-covariance matrix from
#' @param n The number of samples to take (default: 1000)
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `edist` object
#' @export
vcov.edist <- function(object, n = 10000, ...) {
  v <- vcov(ensure_realized(object, n = n))
  # vcov.empirical_dist returns a matrix; drop to scalar for univariate
  drop(v)
}

#' Method for obtaining the mean of an `edist` object.
#' @param x The `edist` object to retrieve the mean from
#' @param n The number of samples to take (default: 1000)
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `edist` object
#' @export
mean.edist <- function(x, n = 10000, ...) {
  mean(ensure_realized(x, n = n))
}

#' Format method for `edist` objects.
#' @param x The object to format
#' @param ... Additional arguments (not used)
#' @return A character string
#' @export
format.edist <- function(x, ...) {
  sprintf("Expression distribution: %s", deparse(x$e))
}

#' Print method for `edist` objects.
#' @param x The object to print
#' @param ... Additional arguments to pass (not used)
#' @export
print.edist <- function(x, ...) {
  cat(format(x), "\n")
  for (i in seq_along(x$vars)) {
    cat("  ", names(x$vars)[i], " ~ ", format(x$vars[[i]]), "\n")
  }
  invisible(x)
}

#' Method for obtaining the sampler of an `edist` object.
#' @param x The `edist` object to obtain the sampler of.
#' @param ... Additional arguments to pass into each of the `sampler`
#'            function generators.
#' @return A function that takes a number of samples `n`, `...`
#'         which is passed into the expression `x$e` and returns
#'         the result of applying the expression `x$e` to the
#'         sampled values.
#' @export
sampler.edist <- function(x, ...) {

  # we use `sampler` to obtain the sampler
  # for each `dist` object in `x$vars`
  # so, here is what `x$vars` looks like:
  # `x$vars <- list("x" = dist1, "y" = dist2)`
  # and here is an example of what `x$e` looks like:
  # `x$e <- expression(x + y)`
  # so, we want to sample from `dist1` and `dist2`
  # and then evaluate `x + y` against the sampled values
  samplers <- lapply(x$vars, sampler, ...)
  cnames <- names(x$vars)

  # todo: accept a param and override the params in the e$vars
  # distributions. each distribution should have a nparam method
  # we can call (including `edist`!), so use it
  function(n = 1, ...) {
    # we sample from each `dist` object
    samples <- sapply(samplers, function(sampler) sampler(n, ...))
    # When n=1, sapply returns a named vector instead of a matrix
    if (!is.matrix(samples)) samples <- matrix(samples, nrow = 1)
    colnames(samples) <- cnames
    esamples <- list()
    for (i in seq_len(nrow(samples))) {
      es <- base::eval(x$e, envir = as.list(samples[i,]))
      esamples[[i]] <- es
    }
    # convert esamples to vector or matrix, depending on
    # dimensionality of the first element in the list
    if (is.matrix(esamples[[1]])) {
      esamples <- do.call(rbind, esamples)
    } else {
      esamples <- unlist(esamples)
    }
    esamples
  }
}



#' Method for adding `dist` objects, or shifting a distribution by a scalar.
#'
#' Creates an expression distribution and automatically simplifies to
#' closed form when possible (e.g., normal + normal = normal,
#' normal + scalar = normal with shifted mean).
#'
#' @param x A `dist` object or numeric scalar
#' @param y A `dist` object or numeric scalar
#' @return A simplified distribution or `edist` if no closed form exists
#' @export
`+.dist` <- function(x, y) {
  if (is.numeric(x) && length(x) == 1) {
    # scalar + dist
    return(simplify(edist(substitute(c + x, list(c = x)), list(x = y))))
  }
  if (is.numeric(y) && length(y) == 1) {
    # dist + scalar
    return(simplify(edist(substitute(x + c, list(c = y)), list(x = x))))
  }
  simplify(edist(quote(x + y), list(x = x, y = y)))
}

#' Method for negation or subtraction of `dist` objects.
#'
#' Unary: returns negated distribution (e.g., -N(mu, var) = N(-mu, var))
#' Binary: creates expression distribution and simplifies to closed form
#' when possible (e.g., normal - normal = normal, normal - scalar = normal).
#'
#' @param x A `dist` object or numeric scalar
#' @param y A `dist` object or numeric scalar (optional for unary negation)
#' @return A simplified distribution or `edist` if no closed form exists
#' @export
`-.dist` <- function(x, y) {
  if (missing(y)) {
    # Unary negation
    if (is_normal(x)) {
      return(normal(mu = -mean(x), var = vcov(x)))
    }
    if (is_uniform_dist(x)) {
      return(uniform_dist(min = -x$max, max = -x$min))
    }
    # For other distributions, create edist with negation
    return(simplify(edist(quote(-x), list(x = x))))
  }
  if (is.numeric(y) && length(y) == 1) {
    # dist - scalar
    return(simplify(edist(substitute(x - c, list(c = y)), list(x = x))))
  }
  if (is.numeric(x) && length(x) == 1) {
    # scalar - dist
    return(simplify(edist(substitute(c - x, list(c = x)), list(x = y))))
  }
  # Binary subtraction
  simplify(edist(quote(x - y), list(x = x, y = y)))
}


#' Multiplication of distribution objects.
#'
#' Handles scalar * dist, dist * scalar, and dist * dist.
#' @param x first operand
#' @param y second operand
#' @return A simplified distribution or edist
#' @export
`*.dist` <- function(x, y) {
  if (is.numeric(x) && length(x) == 1) {
    # scalar * dist
    return(simplify(edist(substitute(c * x, list(c = x)), list(x = y))))
  }
  if (is.numeric(y) && length(y) == 1) {
    # dist * scalar
    return(simplify(edist(substitute(x * c, list(c = y)), list(x = x))))
  }
  # dist * dist
  simplify(edist(quote(x * y), list(x = x, y = y)))
}

#' Power operator for distribution objects.
#'
#' @param x a dist object (base)
#' @param y a numeric scalar (exponent)
#' @return A simplified distribution or edist
#' @export
`^.dist` <- function(x, y) {
  if (is.numeric(y) && length(y) == 1) {
    return(simplify(edist(substitute(x^c, list(c = y)), list(x = x))))
  }
  simplify(edist(quote(x^y), list(x = x, y = y)))
}

#' Division of distribution objects.
#'
#' Handles dist / scalar (delegates to dist * (1/scalar)),
#' scalar / dist, and dist / dist.
#' @param x first operand
#' @param y second operand
#' @return A simplified distribution or edist
#' @export
`/.dist` <- function(x, y) {
  if (is.numeric(y) && length(y) == 1) {
    # dist / scalar → dist * (1/scalar), reuses scalar multiplication rules
    return(x * (1 / y))
  }
  if (is.numeric(x) && length(x) == 1) {
    # scalar / dist
    return(simplify(edist(substitute(c / x, list(c = x)), list(x = y))))
  }
  # dist / dist
  simplify(edist(quote(x / y), list(x = x, y = y)))
}

#' Math group generic for distribution objects.
#'
#' Handles exp(), log(), sqrt(), abs(), cos(), sin(), etc.
#' @param x a dist object
#' @param ... additional arguments
#' @return A simplified distribution or edist
#' @export
Math.dist <- function(x, ...) {
  op <- .Generic
  expr <- substitute(OP(x), list(OP = as.name(op)))
  simplify(edist(expr, list(x = x)))
}

#' Summary group generic for distribution objects.
#'
#' Handles sum(), prod(), min(), max() of distributions.
#' @param ... dist objects
#' @param na.rm ignored
#' @return A simplified distribution or edist
#' @rdname dist_summary_group
#' @importFrom stats setNames
#' @export
Summary.dist <- function(..., na.rm = FALSE) {
  op <- .Generic
  dists <- list(...)

  if (length(dists) == 1) return(dists[[1]])

  if (op == "sum") return(Reduce(`+`, dists))
  if (op == "prod") return(Reduce(`*`, dists))

  if (op == "min") {
    # Check if all are exponential — min of exp is exp with summed rates
    if (all(vapply(dists, is_exponential, logical(1)))) {
      total_rate <- sum(vapply(dists, function(d) d$rate, numeric(1)))
      return(exponential(total_rate))
    }
    # Otherwise build pmin edist
    n <- length(dists)
    var_names <- paste0("x", seq_len(n))
    vars <- setNames(dists, var_names)
    args <- paste(var_names, collapse = ", ")
    expr <- parse(text = paste0("pmin(", args, ")"))[[1]]
    return(simplify(edist(expr, vars)))
  }

  if (op == "max") {
    n <- length(dists)
    var_names <- paste0("x", seq_len(n))
    vars <- setNames(dists, var_names)
    args <- paste(var_names, collapse = ", ")
    expr <- parse(text = paste0("pmax(", args, ")"))[[1]]
    return(simplify(edist(expr, vars)))
  }

  stop("Unsupported Summary operation: ", op)
}

# ---- edist auto-fallback methods via realize() ----------------------------

#' CDF for expression distributions.
#'
#' Falls back to \code{\link{realize}} to materialize the distribution
#' as an \code{\link{empirical_dist}}, then delegates to
#' \code{\link{cdf.empirical_dist}}.
#'
#' @param x An \code{edist} object.
#' @param ... Additional arguments forwarded to \code{cdf.empirical_dist}.
#' @return A function computing the empirical CDF.
#' @export
cdf.edist <- function(x, ...) {
  cdf(ensure_realized(x), ...)
}

#' Density for expression distributions.
#'
#' Falls back to \code{\link{realize}} and delegates to
#' \code{\link{density.empirical_dist}}.
#'
#' @param x An \code{edist} object.
#' @param ... Additional arguments forwarded to \code{density.empirical_dist}.
#' @return A function computing the empirical density (PMF).
#' @export
density.edist <- function(x, ...) {
  density(ensure_realized(x), ...)
}

#' Support for expression distributions.
#'
#' Falls back to \code{\link{realize}} and delegates to
#' \code{\link{sup.empirical_dist}}.
#'
#' @param x An \code{edist} object.
#' @return A \code{finite_set} support object.
#' @export
sup.edist <- function(x) {
  sup(ensure_realized(x))
}

#' Conditional distribution for expression distributions.
#'
#' Falls back to \code{\link{realize}} and delegates to
#' \code{\link{conditional.empirical_dist}}.
#'
#' @param x An \code{edist} object.
#' @param P Predicate function to condition on.
#' @param ... Additional arguments forwarded to the predicate \code{P}.
#' @return A conditional \code{empirical_dist}.
#' @export
conditional.edist <- function(x, P, ...) {
  conditional(ensure_realized(x), P, ...)
}

#' Map function over expression distribution.
#'
#' Falls back to \code{\link{realize}} and delegates to
#' \code{\link{rmap.empirical_dist}}.
#'
#' @param x An \code{edist} object.
#' @param g Function to apply to each observation.
#' @param ... Additional arguments forwarded to \code{g}.
#' @return A transformed \code{empirical_dist}.
#' @export
rmap.edist <- function(x, g, ...) {
  rmap(ensure_realized(x), g, ...)
}

#' Inverse CDF (quantile function) for expression distributions.
#'
#' Falls back to \code{\link{realize}} and delegates to
#' \code{\link{inv_cdf.empirical_dist}}.
#'
#' @param x An \code{edist} object.
#' @param ... Additional arguments forwarded to
#'   \code{inv_cdf.empirical_dist}.
#' @return A function computing the empirical quantile function.
#' @export
inv_cdf.edist <- function(x, ...) {
  inv_cdf(ensure_realized(x), ...)
}