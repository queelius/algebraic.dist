#' Generic method for simplifying distributions.
#'
#' @param x The distribution to simplify
#' @param ... Additional arguments to pass
#' @return The simplified distribution
#' @export
simplify <- function(x, ...) UseMethod("simplify")

#' Default Method for simplifying a `dist` object. Just returns the object.
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
#' when mathematical identities apply. Supported rules include:
#'
#' Single-variable:
#' - c * Normal(mu, v) -> Normal(c*mu, c^2*v)
#' - c * Gamma(a, r) -> Gamma(a, r/c) for c > 0
#' - c * Exponential(r) -> Gamma(1, r/c) for c > 0
#' - Normal(mu, v) + c -> Normal(mu + c, v)
#' - Normal(mu, v) - c -> Normal(mu - c, v)
#' - Normal(0, 1) ^ 2 -> ChiSquared(1)
#' - exp(Normal(mu, v)) -> LogNormal(mu, sqrt(v))
#' - log(LogNormal(ml, sl)) -> Normal(ml, sl^2)
#'
#' Two-variable:
#' - Normal + Normal -> Normal
#' - Normal - Normal -> Normal
#' - Gamma + Gamma (same rate) -> Gamma
#' - Exponential + Exponential (same rate) -> Gamma(2, rate)
#' - Gamma + Exponential (same rate) -> Gamma(a+1, rate)
#' - ChiSquared + ChiSquared -> ChiSquared
#' - Poisson + Poisson -> Poisson
#'
#' @param x The `edist` object to simplify
#' @param ... Additional arguments to pass (not used)
#' @return The simplified distribution, or unchanged `edist` if no rule applies
#' @export
simplify.edist <- function(x, ...) {
  expr <- x$e
  vars <- x$vars

  op <- if (is.call(expr)) as.character(expr[[1]]) else NULL

  # ---- Single-variable rules (length(vars) == 1) ----
  if (length(vars) == 1) {
    d <- vars[[1]]

    # Scalar multiplication: c * dist or dist * c
    if (identical(op, "*")) {
      # Extract scalar from expression
      lhs <- expr[[2]]
      rhs <- expr[[3]]
      scalar <- if (is.numeric(lhs)) lhs else if (is.numeric(rhs)) rhs else NULL

      if (!is.null(scalar)) {
        # c * Normal(mu, v) -> Normal(c*mu, c^2*v)
        if (is_normal(d)) {
          return(normal(mu = scalar * mean(d), var = scalar^2 * vcov(d)))
        }
        # c * Gamma(a, r) -> Gamma(a, r/c) when c > 0
        if (scalar > 0 && is_gamma_dist(d)) {
          return(gamma_dist(shape = d$shape, rate = d$rate / scalar))
        }
        # c * Exponential(r) -> Gamma(1, r/c) when c > 0
        if (scalar > 0 && is_exponential(d)) {
          return(gamma_dist(shape = 1, rate = d$rate / scalar))
        }
      }
    }

    # Location shift: dist + c or c + dist or dist - c
    if (identical(op, "+") || identical(op, "-")) {
      lhs <- expr[[2]]
      rhs <- expr[[3]]

      if (identical(op, "+")) {
        scalar <- if (is.numeric(lhs)) lhs else if (is.numeric(rhs)) rhs else NULL
      } else {
        # For subtraction: x - c, scalar is rhs and negative
        scalar <- if (is.numeric(rhs)) -rhs else NULL
      }

      if (!is.null(scalar) && is_normal(d)) {
        return(normal(mu = mean(d) + scalar, var = vcov(d)))
      }
    }

    # Power: dist ^ c
    if (identical(op, "^")) {
      exponent <- expr[[3]]
      if (is.numeric(exponent) && exponent == 2) {
        # Normal(0,1)^2 -> chi_squared(1)
        if (is_normal(d) && mean(d) == 0 && vcov(d) == 1) {
          return(chi_squared(1))
        }
      }
    }

    # Math functions: exp, log
    if (identical(op, "exp")) {
      # exp(Normal(mu, v)) -> lognormal(mu, sqrt(v))
      if (is_normal(d)) {
        return(lognormal(meanlog = mean(d), sdlog = sqrt(vcov(d))))
      }
    }

    if (identical(op, "log")) {
      # log(LogNormal(ml, sl)) -> Normal(ml, sl^2)
      if (is_lognormal(d)) {
        return(normal(mu = d$meanlog, var = d$sdlog^2))
      }
    }
  }

  # ---- Two-variable rules (length(vars) == 2) ----
  if (length(vars) == 2) {
    d1 <- vars[[1]]
    d2 <- vars[[2]]

    if (identical(op, "+")) {
      # Normal + Normal -> Normal
      if (is_normal(d1) && is_normal(d2)) {
        return(normal(mu = mean(d1) + mean(d2), var = vcov(d1) + vcov(d2)))
      }

      # Gamma + Gamma (same rate) -> Gamma
      if (is_gamma_dist(d1) && is_gamma_dist(d2) && d1$rate == d2$rate) {
        return(gamma_dist(shape = d1$shape + d2$shape, rate = d1$rate))
      }

      # Exponential + Exponential (same rate) -> Gamma(2, rate)
      if (is_exponential(d1) && is_exponential(d2) && d1$rate == d2$rate) {
        return(gamma_dist(shape = 2, rate = d1$rate))
      }

      # Gamma + Exponential (same rate) -> Gamma(a+1, rate)
      if (is_gamma_dist(d1) && is_exponential(d2) && d1$rate == d2$rate) {
        return(gamma_dist(shape = d1$shape + 1, rate = d1$rate))
      }
      if (is_exponential(d1) && is_gamma_dist(d2) && d1$rate == d2$rate) {
        return(gamma_dist(shape = d2$shape + 1, rate = d2$rate))
      }

      # ChiSq + ChiSq -> ChiSq
      if (is_chi_squared(d1) && is_chi_squared(d2)) {
        return(chi_squared(df = d1$df + d2$df))
      }

      # Poisson + Poisson -> Poisson
      if (is_poisson_dist(d1) && is_poisson_dist(d2)) {
        return(poisson_dist(lambda = d1$lambda + d2$lambda))
      }
    }

    if (identical(op, "-")) {
      # Normal - Normal -> Normal
      if (is_normal(d1) && is_normal(d2)) {
        return(normal(mu = mean(d1) - mean(d2), var = vcov(d1) + vcov(d2)))
      }
    }
  }

  # No rule matched
  x
}
