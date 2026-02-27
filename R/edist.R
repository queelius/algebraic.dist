#' Takes an expression `e` and a list `vars` and returns a
#' lazy `edist` (expression distribution object), that is a subclass
#' of `dist` that can be used in place of a `dist` object.
#' 
#' @param e the expression to evaluate against the arguments.
#' @param vars the list of distributions (with variable names)
#'             to evaluate the expression `e` against.
#' @return An `edist` object.
#' @examples
#' x <- normal(0, 1)
#' y <- normal(2, 3)
#' e <- edist(quote(x + y), list(x = x, y = y))
#' e
#' @export
edist <- function(e, vars) {
  primary_classes <- vapply(vars, function(v) class(v)[1], character(1))
  expr_str <- gsub(" ", "_", paste0(e, "_", paste0(primary_classes, collapse = "_")))
  structure(list(e = e, vars = vars, .cache = new.env(parent = emptyenv())),
            class = c(expr_str, "edist", "dist"))
}

#' Function to determine whether an object `x` is an `edist` object.
#' @param x The object to test
#' @return Logical; \code{TRUE} if \code{x} is an \code{edist}.
#' @examples
#' is_edist(normal(0, 1) * exponential(1))  # TRUE
#' is_edist(normal(0, 1))                   # FALSE
#' @export
is_edist <- function(x) {
  inherits(x, "edist")
}

#' Method for obtaining the parameters of an `edist` object.
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @examples
#' z <- normal(0, 1) * exponential(2)
#' params(z)
#' @export
params.edist <- function(x) {
  do.call(c, lapply(x$vars, params))
}

#' Method for obtaining the variance-covariance matrix (or scalar)
#' 
#' @param object The `edist` object to retrieve the variance-covariance matrix from
#' @param n The number of samples to take (default: 10000)
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `edist` object
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(2)
#' vcov(z)
#' }
#' @export
vcov.edist <- function(object, n = 10000, ...) {
  v <- vcov(ensure_realized(object, n = n))
  # vcov.empirical_dist returns a matrix; drop to scalar for univariate
  drop(v)
}

#' Method for obtaining the mean of an `edist` object.
#' @param x The `edist` object to retrieve the mean from
#' @param n The number of samples to take (default: 10000)
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `edist` object
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(2)
#' mean(z)
#' }
#' @export
mean.edist <- function(x, n = 10000, ...) {
  mean(ensure_realized(x, n = n))
}

#' Format method for `edist` objects.
#' @param x The object to format
#' @param ... Additional arguments (not used)
#' @return A character string
#' @examples
#' z <- normal(0, 1) * exponential(2)
#' format(z)
#' @export
format.edist <- function(x, ...) {
  sprintf("Expression distribution: %s", deparse(x$e))
}

#' Print method for `edist` objects.
#' @param x The object to print
#' @param ... Additional arguments to pass (not used)
#' @examples
#' z <- normal(0, 1) * exponential(2)
#' print(z)
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
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(2)
#' s <- sampler(z)
#' samples <- s(100)
#' head(samples)
#' }
#' @export
sampler.edist <- function(x, ...) {
  samplers <- lapply(x$vars, sampler, ...)
  cnames <- names(x$vars)

  function(n = 1, ...) {
    samples <- sapply(samplers, function(s) s(n, ...))
    if (!is.matrix(samples)) samples <- matrix(samples, nrow = 1)
    colnames(samples) <- cnames

    esamples <- vector("list", nrow(samples))
    for (i in seq_len(nrow(samples))) {
      esamples[[i]] <- base::eval(x$e, envir = as.list(samples[i, ]))
    }

    if (is.matrix(esamples[[1]])) do.call(rbind, esamples) else unlist(esamples)
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
#' @examples
#' # Sum of two normals simplifies to a normal
#' z <- normal(0, 1) + normal(2, 3)
#' z  # Normal(mu = 2, var = 4)
#'
#' # Shift a distribution by a constant
#' normal(0, 1) + 5  # Normal(mu = 5, var = 1)
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
#' @examples
#' # Difference of normals simplifies to a normal
#' z <- normal(5, 2) - normal(1, 3)
#' z  # Normal(mu = 4, var = 5)
#'
#' # Unary negation
#' -normal(3, 1)  # Normal(mu = -3, var = 1)
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
#' @examples
#' # Scalar multiplication simplifies for normal
#' z <- 2 * normal(0, 1)
#' z  # Normal(mu = 0, var = 4)
#'
#' # Product of two distributions yields an edist
#' w <- normal(0, 1) * exponential(1)
#' is_edist(w)  # TRUE
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
#' @examples
#' # Standard normal squared yields chi-squared(1)
#' z <- normal(0, 1)^2
#' z
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
#' @examples
#' # Division by scalar reuses multiplication rule
#' z <- normal(0, 4) / 2
#' z  # Normal(mu = 0, var = 1)
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
#' @examples
#' # exp(Normal) simplifies to LogNormal
#' z <- exp(normal(0, 1))
#' z
#'
#' # sqrt of a distribution (no closed-form rule, remains edist)
#' w <- sqrt(exponential(1))
#' is_edist(w)  # TRUE
#' @export
Math.dist <- function(x, ...) {
  op <- .Generic
  expr <- substitute(OP(x), list(OP = as.name(op)))
  simplify(edist(expr, list(x = x)))
}

# Build an edist that applies a parallel function (pmin, pmax) over distributions.
make_elementwise_edist <- function(fn_name, dists) {
  var_names <- paste0("x", seq_along(dists))
  vars <- setNames(dists, var_names)
  expr <- parse(text = paste0(fn_name, "(", paste(var_names, collapse = ", "), ")"))[[1]]
  simplify(edist(expr, vars))
}

#' Summary group generic for distribution objects.
#'
#' Handles sum(), prod(), min(), max() of distributions.
#' @param ... dist objects
#' @param na.rm ignored
#' @return A simplified distribution or edist
#' @examples
#' # sum() reduces via + operator
#' z <- sum(normal(0, 1), normal(2, 3))
#' z  # Normal(mu = 2, var = 4)
#'
#' # min() of exponentials simplifies
#' w <- min(exponential(1), exponential(2))
#' w  # Exponential(rate = 3)
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
    if (all(vapply(dists, is_exponential, logical(1)))) {
      total_rate <- sum(vapply(dists, function(d) d$rate, numeric(1)))
      return(exponential(total_rate))
    }
    return(make_elementwise_edist("pmin", dists))
  }

  if (op == "max") {
    return(make_elementwise_edist("pmax", dists))
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
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(1)
#' Fz <- cdf(z)
#' Fz(0)
#' }
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
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(1)
#' fz <- density(z)
#' }
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
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(1)
#' sup(z)
#' }
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
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) + exponential(1)
#' z_pos <- conditional(z, function(t) t > 2)
#' mean(z_pos)
#' }
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
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(1)
#' abs_z <- rmap(z, abs)
#' mean(abs_z)
#' }
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
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- normal(0, 1) * exponential(1)
#' qz <- inv_cdf(z)
#' qz(0.5)
#' }
#' @export
inv_cdf.edist <- function(x, ...) {
  inv_cdf(ensure_realized(x), ...)
}