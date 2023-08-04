
#' Method for obtaining the inverse cdf for a `univariate_dist` object.
#'
#' We use Newton's method to solve for the root of `cdf(x)(t) - p`.
#'
#' @param x The object to obtain the inverse cdf of.
#' @param ... Additional arguments to pass.
#' @export
inv_cdf.univariate_dist <- function(x, t0 = NULL, eps = 1e-3, ...) {
  # we assume that sup(x) is a contiguous interval
  if (is.null(t0)) {
      t0 <- sampler(x)(1)
  }
  stopifnot(has(Sx, t0))
  function(p, par = params(x), ...) {
    stopifnot(p >= 0 && p <= 1)

    par <- params(x, par)
    Sx <- sup(x, par, ...)
    Fx <- cdf(x, par, ...)
    fx <- density(x, par, ...)

    t1 <- NULL
    repeat {
      alpha <- 1
      repeat
      {
        t1 <- t0 - alpha * (Fx(t0, ...) - p) / fx(t0, ...)
        print(t1)
        if (has(Sx, t1)) {
          break
        }
        alpha <- alpha / 2
      }
      if (abs(t1 - t0) < eps) {
        break
      }
      t0 <- t1
    }
    t1
  }
}

#' Method for obtaining the inverse cdf of an `mvn` object.
#' @param x The object to obtain the inverse cdf of
#' @param ... Additional arguments to pass (not used)
#' @return A function that computes the inverse cdf of the `mvn` distribution, but
#' the quantile produced is the statistical distance (Mahalonbis distance) that
#' defines the surface of the smallest ellipsoid containing `p` probability mass.
#' @importFrom mvtnorm qmvnorm
#' @export
inv_cdf.mvn <- function(x, ...) {
    function(p, mu = x$mu, sigma = x$sigma, log.p = FALSE, ...) {
        q <- qmvnorm(p, mean = mu, sigma = sigma, tail = "lower.tail", ...)
        if (q$message != "Normal Completion") {
            warning(q$message)
        }
        if (log.p) {
            return(log(q$quantile))
        } else {
            return(q$quantile)
        }
    }
}



# Load required libraries
library(mvtnorm)
library(ggplot2)

# Assume we have 2 correlated systems, with given means and covariance matrix
means <- c(1, 1)
cov_mat <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)

# We want to find the combined performance threshold such that the joint failure probability is 0.05
p <- 0.1
q <- qmvnorm(p = p, mean = means, sigma = cov_mat)

# Generate a grid of points
x <- seq(.5, 1.5, length.out = 100)
y <- seq(.5, 1.5, length.out = 100)
grid <- expand.grid(x = x, y = y)

# For each point in the grid, calculate its Mahalanobis distance from the mean
mahalanobis_distances <- mahalanobis(x = grid, center = means, cov = solve(cov_mat))

# Create a data frame with the grid points and their Mahalanobis distances
grid$dist <- mahalanobis_distances

# Plot the points in the grid, coloring points inside the ellipsoid differently
ggplot(grid, aes(x = x, y = y, color = dist <= q)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(x = "System 1 Performance", y = "System 2 Performance", 
       title = "Performance Ellipsoid for Joint Failure Probability p = 0.05") +
  theme_minimal()









#' This is the expectation of a function `g` with respect to an
#' normal distribution `x`.
#'
#' @param x The disrtibution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass into `g`
#' @importFrom stats integrate
expectation.normal <- function(x, g, ...) {

  stopifnot(is.function(g))
  f <- density(x)
  integrate(function(t) g(t, ...) * f(t), -Inf, Inf)
}




















#' Computes the variance-covariance matrix of a distribution-like object.
#'
#' @param x the distribution-like object to compute the variance-covariance
#' @export
vcov <- function(x) {
    UseMethod("vcov", x)
}

#' Computes the stadnard deviation of a distribution-like object.
#'
#' @param x the distribution-like object to compute the sd of
#' @export
sd <- function(x) {
    UseMethod("sd", x)
}
