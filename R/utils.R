#' Function used for computing expectations given data (e.g., from an MC
#' simulation or bootstrap). it expects a matrix, or something that can be
#' coerced to a matrix (e.g., a data frame). it also expects a function `g` to
#' apply to each row of the data, and returns the expectation of `g` under the
#' empirical distribution of the data. it also returns a confidence interval for
#' the expectation, and the number of samples used to compute the expectation.
#' 
#' example: expectation_data(D, function(x) (x-colMeans(D)) %*% t(x-colMeans(D)))
#' computes the covariance of the data D, except the matrix structure is lost
#' (it's just a vector, which can be coerced back to a matrix if needed).
#' 
#' @importFrom stats qnorm cov var
#' @param data a matrix of data
#' @param g a function to apply to each row of the data
#' @param ... additional arguments to pass to `g`
#' @param compute_stats whether to compute CIs for the expectations
#' @param alpha the confidence level for the confidence interval for each
#'              component of the expectation (if compute_stats is TRUE)
#' @return if compute_stats is TRUE, then a list with the following components:
#'  value - The estimate of the expectation
#'  ci    - The confidence intervals for each component of the expectation
#'  n     - The number of samples
#' otherwise, just the value of the expectation.
#' @export
expectation_data <- function(
  data,
  g = function(x) x,
  ...,
  compute_stats = TRUE,
  alpha = 0.05) {

    if (!is.matrix(data)) {
        data <- as.matrix(data)
    }

    stopifnot(
      is.function(g),
      is.numeric(alpha), alpha > 0, alpha < 1)

    n <- nrow(data)
    data <- apply(data, 1, g, ...)
    if (is.matrix(data)) {
      data <- t(data)
      mu <- colMeans(data)
    } else {
      mu <- mean(data)
    }
    if (!compute_stats) {
      return(mu)
    }

    z.alpha <- qnorm(1 - alpha / 2)
    if (is.matrix(data)) {
        v <- diag(cov(data))
        sd <- sqrt(pmax(v, 0) / n)
        ci <- matrix(NA, nrow = length(mu), ncol = 2)
        ci[, 1] <- mu - z.alpha * sd
        ci[, 2] <- mu + z.alpha * sd
        return(list(value = mu, ci = ci, n = n))
    } else {
        v <- var(data)
        sd <- sqrt(max(v, 0) / n)
        return(list(value = mu,
                    ci = c(mu - z.alpha * sd, mu + z.alpha * sd),
                    n = n))
    }
}
