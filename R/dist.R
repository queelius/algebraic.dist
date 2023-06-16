#' Print method for `summary_dist` objects.
#' @param x The object to print
#' @param ... Additional arguments
#' @export
print.summary_dist <- function(x, ...) {
  cat(x$name, "\n")
  cat("Mean:", x$mean, "\n")
  cat("Variance:", x$vcov, "\n")
}

#' Method for constructing a `summary_dist` object.
#' @param name The name of the distribution
#' @param mean The mean of the distribution
#' @param vcov The variance of the distribution
#' @param ... Additional arguments
#' @return A `summary_dist` object
#' @export
summary_dist <- function(name, mean, vcov, ...) {
  structure(list(
    name = name,
    mean = mean,
    vcov = vcov),
    class = "summary_dist")
}


#' Function to determine whether an object `x` is a `dist` object.
#' @param x The object to test
#' @export
is_dist <- function(x) {
  inherits(x, "dist")
}

#' Expectation operator applied to `object` of type `dist`
#' with respect to a function `g`. That is, `E(g(object))`.
#' 
#' We use the CLT to construct a CI(`alpha`) for the estimate of the
#' expectation. That is, we estimate `E(g(object))` with the sample mean
#' and Var(g(object)) with the sigma^2/n, where sigma^2 is the sample
#' variance of g(object) and n is the number of samples. From these,
#' we construct the CI.
#'
#' @param x `dist` object
#' @param g characteristic function of interest, defaults to identity
#' @param n number of sample points, larger samples reduce variance
#' @param alpha significance level for the confidence interval
#' @param ... additional arguments to pass to `g`
#' @return A list with the following components:
#'   value - The estimate of the expectation
#'   ci    - The confidence interval
#'   n     - The number of samples
#' @importFrom stats qnorm cov var
#' @export
expectation.dist <- function(
    x,
    g = function(t) t,
    ...,
    n = 10000L,
    alpha = 0.05) {

    stopifnot(is.function(g), is.numeric(n), n > 0,
              is.numeric(alpha), alpha > 0, alpha < 1)

    data <- sampler(x)(n)
    data <- if (is.matrix(data)) {
      apply(data, 1, g, ...)
    } else {
      g(data)
    }

    z.alpha <- qnorm(1 - alpha / 2)
    if (is.matrix(data)) {
        mu <- colMeans(data)
        sd <- sqrt(diag(cov(data)) / n)
        # we construct a confidence interval for
        # each component of the expectation
        ci <- matrix(NA, nrow = length(mu), ncol = 2)
        ci[, 1] <- mu - z.alpha * sd
        ci[, 2] <- mu + z.alpha * sd
        return(list(value = mu,
                    ci = ci,
                    n = n))
    } else {
        mu <- mean(data)
        sd <- sqrt(var(data) / n)
        return(list(value = mu,
                    ci = c(mu - z.alpha * sd, mu + z.alpha * sd),
                    n = n))
    }
  }
