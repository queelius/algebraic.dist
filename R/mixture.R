#' Construct a mixture distribution.
#'
#' Creates an S3 object representing a finite mixture distribution.
#' The density is \eqn{f(x) = \sum_{k=1}^{K} w_k f_k(x)} where
#' \eqn{f_k} are the component densities and \eqn{w_k} are the mixing
#' weights.
#'
#' The class hierarchy is determined by the components: if all components
#' are univariate (or multivariate, continuous, discrete), the mixture
#' inherits those classes as well.
#'
#' @param components A non-empty list of \code{dist} objects.
#' @param weights A numeric vector of non-negative mixing weights that
#'   sum to 1 (within tolerance \code{1e-10}).  Must have the same
#'   length as \code{components}.
#' @return A \code{mixture} object with appropriate class hierarchy.
#' @export
mixture <- function(components, weights) {
  if (!is.list(components) || length(components) == 0)
    stop("'components' must be a non-empty list of distributions")
  if (!all(sapply(components, inherits, "dist")))
    stop("all components must be 'dist' objects")
  if (!is.numeric(weights) || length(weights) != length(components))
    stop("'weights' must be a numeric vector of length ", length(components))
  if (any(weights < 0))
    stop("'weights' must be non-negative")
  if (abs(sum(weights) - 1) > 1e-10)
    stop("'weights' must sum to 1, got: ", sum(weights))

  # Determine class based on components
  classes <- "mixture"
  if (all(sapply(components, inherits, "univariate_dist"))) {
    classes <- c(classes, "univariate_dist")
  } else if (all(sapply(components, inherits, "multivariate_dist"))) {
    classes <- c(classes, "multivariate_dist")
  }
  if (all(sapply(components, inherits, "continuous_dist"))) {
    classes <- c(classes, "continuous_dist")
  } else if (all(sapply(components, inherits, "discrete_dist"))) {
    classes <- c(classes, "discrete_dist")
  }
  classes <- c(classes, "dist")

  structure(list(components = components, weights = weights),
            class = classes)
}

#' Test whether an object is a \code{mixture} distribution.
#'
#' @param x The object to test.
#' @return \code{TRUE} if \code{x} inherits from \code{"mixture"},
#'   \code{FALSE} otherwise.
#' @export
is_mixture <- function(x) inherits(x, "mixture")

#' Mean of a mixture distribution.
#'
#' The mean of a mixture is the weighted sum of the component means:
#' \eqn{E[X] = \sum_k w_k \mu_k}.
#'
#' @param x A \code{mixture} object.
#' @param ... Additional arguments (not used).
#' @return The mean of the mixture distribution.
#' @export
mean.mixture <- function(x, ...) {
  weighted_means <- mapply(function(comp, w) w * mean(comp),
                           x$components, x$weights, SIMPLIFY = FALSE)
  Reduce(`+`, weighted_means)
}

#' Variance of a mixture distribution.
#'
#' Uses the law of total variance:
#' \eqn{Var(X) = E[Var(X|K)] + Var(E[X|K])}.
#'
#' @param object A \code{mixture} object.
#' @param ... Additional arguments (not used).
#' @return The variance (scalar for univariate mixtures).
#' @export
vcov.mixture <- function(object, ...) {
  d <- dim(object$components[[1]])
  if (d == 1) {
    # Univariate law of total variance
    comp_means <- sapply(object$components, mean)
    comp_vars <- sapply(object$components, vcov)
    overall_mean <- sum(object$weights * comp_means)
    within_var <- sum(object$weights * comp_vars)
    between_var <- sum(object$weights * (comp_means - overall_mean)^2)
    within_var + between_var
  } else {
    # Multivariate law of total variance
    overall_mean <- mean(object)
    within_cov <- matrix(0, d, d)
    between_cov <- matrix(0, d, d)
    for (k in seq_along(object$components)) {
      w <- object$weights[k]
      within_cov <- within_cov + w * vcov(object$components[[k]])
      dev <- mean(object$components[[k]]) - overall_mean
      between_cov <- between_cov + w * outer(dev, dev)
    }
    within_cov + between_cov
  }
}

#' Probability density function for a mixture distribution.
#'
#' Returns a function that evaluates the mixture density at given points.
#' The mixture density is \eqn{f(x) = \sum_k w_k f_k(x)}.
#'
#' @param x A \code{mixture} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(t, log = FALSE, ...)} returning the
#'   density (or log-density) at \code{t}.
#' @importFrom stats density
#' @export
density.mixture <- function(x, ...) {
  comp_densities <- lapply(x$components, density)
  function(t, log = FALSE, ...) {
    vals <- sapply(seq_along(x$components), function(i) {
      x$weights[i] * comp_densities[[i]](t, ...)
    })
    result <- if (is.matrix(vals)) rowSums(vals) else sum(vals)
    if (log) log(result) else result
  }
}

#' Cumulative distribution function for a mixture distribution.
#'
#' Returns a function that evaluates the mixture CDF at given points.
#' The mixture CDF is \eqn{F(x) = \sum_k w_k F_k(x)}.
#'
#' @param x A \code{mixture} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(q, ...)} returning the CDF at \code{q}.
#' @export
cdf.mixture <- function(x, ...) {
  comp_cdfs <- lapply(x$components, cdf)
  function(q, ...) {
    vals <- sapply(seq_along(x$components), function(i) {
      x$weights[i] * comp_cdfs[[i]](q, ...)
    })
    if (is.matrix(vals)) rowSums(vals) else sum(vals)
  }
}

#' Sampler for a mixture distribution.
#'
#' Returns a function that draws samples from the mixture by first
#' selecting a component according to the mixing weights, then
#' sampling from the selected component.
#'
#' @param x A \code{mixture} object.
#' @param ... Additional arguments (not used).
#' @return A function \code{function(n = 1, ...)} returning a numeric
#'   vector of length \code{n}.
#' @export
sampler.mixture <- function(x, ...) {
  comp_samplers <- lapply(x$components, sampler)
  d <- dim(x$components[[1]])
  is_mv <- d > 1
  function(n = 1, ...) {
    # Draw component indices
    indices <- sample(length(x$components), n, replace = TRUE,
                      prob = x$weights)
    if (is_mv) {
      result <- matrix(NA_real_, nrow = n, ncol = d)
      for (k in seq_along(x$components)) {
        mask <- indices == k
        if (any(mask)) {
          result[mask, ] <- comp_samplers[[k]](sum(mask), ...)
        }
      }
    } else {
      result <- numeric(n)
      for (k in seq_along(x$components)) {
        mask <- indices == k
        if (any(mask)) {
          result[mask] <- comp_samplers[[k]](sum(mask), ...)
        }
      }
    }
    result
  }
}

#' Retrieve the parameters of a \code{mixture} object.
#'
#' Returns a named numeric vector containing all component parameters
#' (flattened) followed by the mixing weights.
#'
#' @param x A \code{mixture} object.
#' @return A named numeric vector.
#' @export
params.mixture <- function(x) {
  comp_params <- lapply(x$components, params)
  all_params <- unlist(comp_params)
  c(all_params, setNames(x$weights,
    paste0("weight", seq_along(x$weights))))
}

#' Number of parameters for a \code{mixture} distribution.
#'
#' The total number of parameters is the sum of component parameters
#' plus the number of mixing weights.
#'
#' @param x A \code{mixture} object.
#' @return An integer count of parameters.
#' @export
nparams.mixture <- function(x) {
  sum(sapply(x$components, function(comp) length(params(comp)))) +
    length(x$weights)
}

#' Dimension of a mixture distribution.
#'
#' Returns the dimension of the first component (all components are
#' assumed to have the same dimension).
#'
#' @param x A \code{mixture} object.
#' @return The dimension of the distribution.
#' @export
dim.mixture <- function(x) {
  dim(x$components[[1]])
}

#' Support of a mixture distribution.
#'
#' Returns an \code{\link{interval}} spanning the widest range of all
#' component supports (from the smallest infimum to the largest supremum).
#'
#' @param x A \code{mixture} object.
#' @return An \code{interval} object.
#' @export
sup.mixture <- function(x) {
  sups <- lapply(x$components, sup)
  lo <- min(sapply(sups, infimum))
  hi <- max(sapply(sups, supremum))
  interval$new(lower = lo, upper = hi)
}

#' Format a \code{mixture} object as a character string.
#'
#' @param x A \code{mixture} object.
#' @param ... Additional arguments (not used).
#' @return A character string describing the mixture.
#' @export
format.mixture <- function(x, ...) {
  k <- length(x$components)
  sprintf("Mixture distribution (%d components)", k)
}

#' Print a \code{mixture} object.
#'
#' @param x A \code{mixture} object.
#' @param ... Additional arguments (not used).
#' @return \code{x}, invisibly.
#' @export
print.mixture <- function(x, ...) {
  cat(format(x), "\n")
  for (i in seq_along(x$components)) {
    cat(sprintf("  [w=%.3f] %s\n", x$weights[i], format(x$components[[i]])))
  }
  invisible(x)
}


#' Marginal distribution of a mixture.
#'
#' The marginal of a mixture is itself a mixture of the component marginals
#' with the same mixing weights:
#' \eqn{p(x_I) = \sum_k w_k p_k(x_I)}.
#'
#' Requires all components to support \code{\link{marginal}}.
#'
#' @param x A \code{mixture} object.
#' @param indices Integer vector of variable indices to keep.
#' @return A \code{mixture} object with marginalized components.
#' @export
marginal.mixture <- function(x, indices) {
  comp_marginals <- lapply(x$components, marginal, indices = indices)
  mixture(comp_marginals, x$weights)
}


#' Conditional distribution of a mixture.
#'
#' For a mixture of distributions that support closed-form conditioning
#' (e.g. MVN), uses Bayes' rule to update the mixing weights:
#' \deqn{w_k' \propto w_k f_k(x_{given})}
#' where \eqn{f_k} is the marginal density of component \eqn{k} at the
#' observed values. The component conditionals are computed via
#' \code{conditional(component_k, given_indices = ..., given_values = ...)}.
#'
#' Falls back to MC realization if \code{P} is provided or if any
#' component does not support \code{given_indices}/\code{given_values}.
#'
#' @param x A \code{mixture} object.
#' @param P Optional predicate function for MC fallback.
#' @param ... Additional arguments.
#' @param given_indices Integer vector of observed variable indices.
#' @param given_values Numeric vector of observed values.
#' @return A \code{mixture} or \code{empirical_dist} object.
#' @export
conditional.mixture <- function(x, P = NULL, ...,
                                given_indices = NULL, given_values = NULL) {
  # Closed-form path for mixture-of-MVN
  if (!is.null(given_indices) && !is.null(given_values)) {
    K <- length(x$components)

    # Compute marginal densities at observed values for weight update
    log_weights <- numeric(K)
    for (k in seq_len(K)) {
      comp <- x$components[[k]]
      # Get marginal density over given_indices
      marg <- marginal(comp, given_indices)
      dens_fn <- density(marg)
      log_weights[k] <- log(x$weights[k]) + log(dens_fn(given_values))
    }

    # Normalize weights (log-sum-exp for numerical stability)
    max_lw <- max(log_weights)
    new_weights <- exp(log_weights - max_lw)
    new_weights <- new_weights / sum(new_weights)

    # Compute component conditionals
    new_components <- lapply(x$components, function(comp) {
      conditional(comp, given_indices = given_indices,
                  given_values = given_values)
    })

    return(mixture(new_components, new_weights))
  }

  # Predicate-based MC fallback
  if (!is.null(P)) {
    return(conditional(ensure_realized(x), P, ...))
  }

  stop("must provide either 'P' or both 'given_indices' and 'given_values'")
}
