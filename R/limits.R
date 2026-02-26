#' Central Limit Theorem Limiting Distribution
#'
#' Returns the limiting distribution of the standardized sample mean
#' \eqn{\sqrt{n}(\bar{X}_n - \mu)} under the Central Limit Theorem.
#' For a univariate distribution with variance \eqn{\sigma^2}, this is
#' \eqn{N(0, \sigma^2)}. For a multivariate distribution with covariance
#' matrix \eqn{\Sigma}, this is \eqn{MVN(0, \Sigma)}.
#'
#' @param base_dist A \code{dist} object representing the base distribution.
#' @return A \code{normal} or \code{mvn} distribution representing the
#'   CLT limiting distribution.
#' @importFrom stats vcov
#' @export
clt <- function(base_dist) {
  if (!inherits(base_dist, "dist"))
    stop("'base_dist' must be a 'dist' object, got class: ",
         paste(class(base_dist), collapse = ", "))

  d <- dim(base_dist)
  v <- vcov(base_dist)

  if (is.null(d) || d == 1) {
    normal(mu = 0, var = v)
  } else {
    mvn(mu = rep(0, d), sigma = v)
  }
}

#' Law of Large Numbers Limiting Distribution
#'
#' Returns the degenerate limiting distribution of the sample mean
#' \eqn{\bar{X}_n} under the Law of Large Numbers. The limit is a
#' point mass at the population mean (represented as a normal or mvn
#' with zero variance).
#'
#' @param base_dist A \code{dist} object representing the base distribution.
#' @return A \code{normal} or \code{mvn} distribution with zero variance,
#'   representing the degenerate distribution at the mean.
#' @importFrom stats vcov
#' @export
lln <- function(base_dist) {
  if (!inherits(base_dist, "dist"))
    stop("'base_dist' must be a 'dist' object, got class: ",
         paste(class(base_dist), collapse = ", "))

  d <- dim(base_dist)
  mu <- mean(base_dist)

  if (is.null(d) || d == 1) {
    normal(mu = mu, var = 0)
  } else {
    mvn(mu = mu, sigma = matrix(0, d, d))
  }
}

#' Delta Method CLT Limiting Distribution
#'
#' Returns the limiting distribution of \eqn{\sqrt{n}(g(\bar{X}_n) - g(\mu))}
#' under the Delta Method. For a univariate distribution, this is
#' \eqn{N(0, g'(\mu)^2 \sigma^2)}. For a multivariate distribution with
#' Jacobian \eqn{J = Dg(\mu)}, this is \eqn{MVN(0, J \Sigma J^T)}.
#'
#' @param base_dist A \code{dist} object representing the base distribution.
#' @param g The function to apply to the sample mean.
#' @param dg The derivative (univariate) or Jacobian function (multivariate)
#'   of \code{g}. For univariate distributions, \code{dg(x)} should return
#'   a scalar. For multivariate distributions, \code{dg(x)} should return
#'   a matrix (the Jacobian).
#' @return A \code{normal} or \code{mvn} distribution representing the
#'   Delta Method limiting distribution.
#' @importFrom stats vcov
#' @export
delta_clt <- function(base_dist, g, dg) {
  if (!inherits(base_dist, "dist"))
    stop("'base_dist' must be a 'dist' object, got class: ",
         paste(class(base_dist), collapse = ", "))

  d <- dim(base_dist)
  mu <- mean(base_dist)
  v <- vcov(base_dist)

  if (is.null(d) || d == 1) {
    dg_val <- dg(mu)
    normal(mu = 0, var = dg_val^2 * v)
  } else {
    J <- dg(mu)
    if (!is.matrix(J))
      stop("'dg' must return a matrix (Jacobian) for multivariate distributions")
    sigma_out <- J %*% v %*% t(J)
    mvn(mu = rep(0, nrow(J)), sigma = sigma_out)
  }
}

#' Moment-Matching Normal Approximation
#'
#' Constructs a normal (or multivariate normal) distribution that matches
#' the mean and variance-covariance of the input distribution. This is
#' useful as a quick Gaussian approximation for any distribution whose
#' first two moments are available.
#'
#' @param x A \code{dist} object to approximate.
#' @return A \code{normal} distribution (for univariate inputs) or an
#'   \code{mvn} distribution (for multivariate inputs) with the same
#'   mean and variance-covariance as \code{x}.
#' @importFrom stats vcov
#' @export
normal_approx <- function(x) {
  if (!inherits(x, "dist"))
    stop("'x' must be a 'dist' object, got class: ",
         paste(class(x), collapse = ", "))

  d <- dim(x)
  mu <- mean(x)
  v <- vcov(x)

  if (is.null(d) || d == 1) {
    normal(mu = mu, var = v)
  } else {
    mvn(mu = mu, sigma = v)
  }
}
