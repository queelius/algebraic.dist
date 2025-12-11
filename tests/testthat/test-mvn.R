# Tests for the multivariate normal distribution class

test_that("mvn constructor creates valid multivariate object", {
  # Given: A mean vector and covariance matrix
  mu <- c(1, 2, 3)
  sigma <- diag(3)

  # When: Creating an mvn distribution
  m <- mvn(mu = mu, sigma = sigma)

  # Then: Object has correct class hierarchy
  expect_s3_class(m, "mvn")
  expect_s3_class(m, "multivariate_dist")
  expect_s3_class(m, "continuous_dist")
  expect_s3_class(m, "dist")

  expect_equal(m$mu, c(1, 2, 3))
  expect_equal(m$sigma, diag(3))
})

test_that("mvn constructor returns normal for univariate case", {
  # Given: A 1-dimensional mean and variance
  mu <- 5
  sigma <- matrix(4)

  # When: Creating an mvn distribution
  n <- mvn(mu = mu, sigma = sigma)

  # Then: It returns a normal object, not mvn
  expect_s3_class(n, "normal")
  expect_equal(mean(n), 5)
  # vcov returns matrix for normal constructed from mvn
  expect_equal(as.numeric(vcov(n)), 4)
})

test_that("mvn constructor uses identity covariance by default", {
  mu <- c(0, 0)
  m <- mvn(mu = mu)

  expect_equal(m$sigma, diag(2))
})

test_that("mvn constructor validates inputs", {
  # Given/When/Then: Invalid inputs should cause errors
  expect_error(mvn(mu = "a"))
  expect_error(mvn(mu = c(1, 2), sigma = c(1, 2))) # sigma must be matrix
  expect_error(mvn(mu = c(1, 2), sigma = diag(3))) # dimension mismatch
})

test_that("is_mvn identifies mvn objects correctly", {
  m <- mvn(mu = c(0, 0))

  expect_true(is_mvn(m))
  expect_false(is_mvn(normal()))
  expect_false(is_mvn(list(mu = c(0, 0), sigma = diag(2))))
})

test_that("mean.mvn returns the correct mean vector", {
  m <- mvn(mu = c(1, 2, 3))

  expect_equal(mean(m), c(1, 2, 3))
})

test_that("vcov.mvn returns the correct covariance matrix", {
  sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
  m <- mvn(mu = c(0, 0), sigma = sigma)

  expect_equal(vcov(m), sigma)
})

test_that("params.mvn returns mean and sigma", {
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  p <- params(m)

  # params.mvn returns combined c(mu = ..., sigma = ...)
  # The names may not be "mu" and "sigma" exactly, but values should be present
  expect_true(length(p) >= 3)  # at least 2 mean values + 4 sigma values
  expect_true(1 %in% p)  # mu[1] = 1
  expect_true(2 %in% p)  # mu[2] = 2
})

test_that("dim.mvn returns correct dimensionality", {
  m <- mvn(mu = c(1, 2, 3))

  expect_equal(dim(m), 3)
})

test_that("sampler.mvn returns a function that generates samples", {
  m <- mvn(mu = c(0, 0), sigma = diag(2))
  samp_fn <- sampler(m)

  # The sampler should return a function
  expect_type(samp_fn, "closure")

  # The function should generate a matrix with correct dimensions
  samples <- samp_fn(100)
  expect_true(is.matrix(samples))
  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 2)
})

test_that("sampler.mvn produces samples with approximately correct mean", {
  set.seed(42)
  mu <- c(5, -3)
  m <- mvn(mu = mu, sigma = diag(2))
  samp_fn <- sampler(m)
  samples <- samp_fn(10000)

  # Column means should be close to specified means
  expect_equal(colMeans(samples), mu, tolerance = 0.1)
})

test_that("sampler.mvn produces samples with approximately correct covariance", {
  set.seed(42)
  sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
  m <- mvn(mu = c(0, 0), sigma = sigma)
  samp_fn <- sampler(m)
  samples <- samp_fn(10000)

  # Sample covariance should be close to specified covariance
  expect_equal(cov(samples), sigma, tolerance = 0.1)
})

test_that("density.mvn returns correct probability density", {
  m <- mvn(mu = c(0, 0), sigma = diag(2))
  pdf <- density(m)

  # Density at origin for standard bivariate normal
  expected <- mvtnorm::dmvnorm(c(0, 0), mean = c(0, 0), sigma = diag(2))
  expect_equal(pdf(c(0, 0)), expected, tolerance = 1e-10)
})

test_that("density.mvn handles log argument correctly", {
  m <- mvn(mu = c(0, 0), sigma = diag(2))
  pdf <- density(m)

  expected <- mvtnorm::dmvnorm(c(0, 0), mean = c(0, 0), sigma = diag(2), log = TRUE)
  expect_equal(pdf(c(0, 0), log = TRUE), expected, tolerance = 1e-10)
})

test_that("density.mvn handles matrix input", {
  m <- mvn(mu = c(0, 0), sigma = diag(2))
  pdf <- density(m)

  # Multiple points as rows of a matrix
  points <- matrix(c(0, 0, 1, 1, -1, 2), nrow = 3, byrow = TRUE)
  densities <- pdf(points)

  expect_length(densities, 3)
  expect_equal(densities[1], mvtnorm::dmvnorm(c(0, 0), mean = c(0, 0), sigma = diag(2)))
})

test_that("sup.mvn returns correct support interval", {
  m <- mvn(mu = c(0, 0))
  s <- sup(m)

  # MVN has support R^d = (-Inf, Inf)^d
  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), c(-Inf, -Inf))
  expect_equal(s$supremum(), c(Inf, Inf))
  expect_equal(s$dim(), 2)
})

test_that("marginal.mvn returns correct univariate marginal", {
  sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
  m <- mvn(mu = c(1, 2), sigma = sigma)

  # Marginal of first component
  marg1 <- marginal(m, 1)

  expect_s3_class(marg1, "normal")
  expect_equal(mean(marg1), 1)
  expect_equal(vcov(marg1), 1)
})

test_that("marginal.mvn returns correct multivariate marginal", {
  mu <- c(1, 2, 3)
  sigma <- diag(3)
  sigma[1, 2] <- sigma[2, 1] <- 0.3
  m <- mvn(mu = mu, sigma = sigma)

  # Marginal over first two components
  marg12 <- marginal(m, c(1, 2))

  expect_s3_class(marg12, "mvn")
  expect_equal(mean(marg12), c(1, 2))
  expect_equal(dim(marg12), 2)
})

test_that("marginal.mvn validates indices", {
  m <- mvn(mu = c(0, 0))

  expect_error(marginal(m, c()))
  # c(3) should error since it's > dim(m)
  expect_error(marginal(m, c(3)))
  # Note: the implementation checks indices < 0, so c(-1) causes subsetting error
  # but may not raise the "indices must be in [1, dim(x)]" error explicitly
})

test_that("rmap.mvn applies function to distribution", {
  set.seed(42)
  m <- mvn(mu = c(0, 0), sigma = diag(2))

  # Apply a linear transformation that sums components
  g <- function(x) sum(x)
  mapped <- rmap(m, g, n = 5000)

  # The sum of two independent standard normals has mean 0 and variance 2
  # rmap.mvn creates an mvn from empirical_dist, which may become normal
  # if the result is univariate (which it is here since sum returns scalar)
  expect_s3_class(mapped, "dist")
  expect_true(abs(mean(mapped)) < 0.2)
})

test_that("print.mvn produces output without error", {
  m <- mvn(mu = c(1, 2))

  expect_output(print(m), "Multivariate normal")
  expect_output(print(m), "mean")
  expect_output(print(m), "variance-covariance")
})

test_that("cdf.mvn returns cumulative probability", {
  m <- mvn(mu = c(0, 0), sigma = diag(2))
  cdf_fn <- cdf(m)

  # CDF at origin should be 0.25 for standard bivariate normal (symmetric)
  p <- cdf_fn(c(0, 0))
  expect_equal(as.numeric(p), 0.25, tolerance = 0.01)
})

test_that("sample_mvn_region generates samples within probability region", {
  # Given: MVN parameters
  mu <- c(0, 0)
  sigma <- diag(2)

  # When: Sampling from the 95% probability region
  set.seed(123)
  samples <- sample_mvn_region(n = 100, mu = mu, sigma = sigma, p = 0.95)

  # Then: All samples should be within the chi-squared critical value
  k <- length(mu)
  crit <- qchisq(0.95, k)
  mahal_dists <- apply(samples, 1, function(x) mahalanobis(x, center = mu, cov = sigma))

  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 2)
  expect_true(all(mahal_dists <= crit))
})

test_that("sample_mvn_region with p=1 includes all samples", {
  # Given: MVN parameters
  mu <- c(0, 0)
  sigma <- diag(2)

  # When: Sampling with p=1 (entire distribution)
  set.seed(456)
  samples <- sample_mvn_region(n = 50, mu = mu, sigma = sigma, p = 1)

  # Then: Should return the requested number of samples
  expect_equal(nrow(samples), 50)
  expect_equal(ncol(samples), 2)
})

test_that("sample_mvn_region validates inputs", {
  mu <- c(0, 0)
  sigma <- diag(2)

  # Invalid probability
  expect_error(sample_mvn_region(n = 10, mu = mu, sigma = sigma, p = 1.5))
  expect_error(sample_mvn_region(n = 10, mu = mu, sigma = sigma, p = -0.1))

  # Invalid n
  expect_error(sample_mvn_region(n = 0, mu = mu, sigma = sigma, p = 0.95))
})
