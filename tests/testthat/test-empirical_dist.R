# Tests for the empirical distribution class

test_that("empirical_dist constructor creates valid univariate object from vector", {
  # Given: A numeric vector
  data <- c(1, 2, 3, 4, 5)

  # When: Creating an empirical distribution
  e <- empirical_dist(data)

  # Then: Object has correct class hierarchy
  expect_s3_class(e, "empirical_dist")
  expect_s3_class(e, "univariate_dist")
  expect_s3_class(e, "discrete_dist")
  expect_s3_class(e, "dist")
})

test_that("empirical_dist constructor creates valid multivariate object from matrix", {
  # Given: A numeric matrix
  data <- matrix(1:6, nrow = 3, ncol = 2)

  # When: Creating an empirical distribution
  e <- empirical_dist(data)

  # Then: Object has correct class hierarchy
  expect_s3_class(e, "empirical_dist")
  expect_s3_class(e, "multivariate_dist")
  expect_s3_class(e, "discrete_dist")
  expect_s3_class(e, "dist")
})

test_that("empirical_dist constructor validates non-empty data", {
  expect_error(empirical_dist(c()))
  expect_error(empirical_dist(numeric(0)))
})

test_that("is_empirical_dist identifies empirical_dist objects correctly", {
  e <- empirical_dist(1:5)

  expect_true(is_empirical_dist(e))
  expect_false(is_empirical_dist(list(data = 1:5)))
  expect_false(is_empirical_dist(normal()))
})

test_that("dim.empirical_dist returns correct dimensionality", {
  # Univariate
  e1 <- empirical_dist(1:5)
  expect_equal(dim(e1), 1)

  # Multivariate
  e2 <- empirical_dist(matrix(1:10, nrow = 5, ncol = 2))
  expect_equal(dim(e2), 2)
})

test_that("nobs.empirical_dist returns correct number of observations", {
  e <- empirical_dist(1:10)
  expect_equal(nobs(e), 10)

  e2 <- empirical_dist(matrix(1:20, nrow = 5, ncol = 4))
  expect_equal(nobs(e2), 5)
})

test_that("obs.empirical_dist returns the underlying data", {
  data <- c(1, 2, 3, 4, 5)
  e <- empirical_dist(data)

  result <- obs(e)
  expect_equal(as.vector(result), data)
})

test_that("mean.empirical_dist returns correct mean for univariate", {
  data <- c(1, 2, 3, 4, 5)
  e <- empirical_dist(data)

  expect_equal(mean(e), 3)
})

test_that("mean.empirical_dist returns correct column means for multivariate", {
  data <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, ncol = 2)
  e <- empirical_dist(data)

  expect_equal(mean(e), c(2, 20))
})

test_that("vcov.empirical_dist returns correct variance for univariate", {
  data <- c(1, 2, 3, 4, 5)
  e <- empirical_dist(data)

  # vcov returns a matrix even for univariate data (cov() behavior)
  expect_equal(as.numeric(vcov(e)), var(data), tolerance = 1e-10)
})

test_that("vcov.empirical_dist returns correct covariance matrix for multivariate", {
  data <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, ncol = 2)
  e <- empirical_dist(data)

  expect_equal(vcov(e), cov(data), tolerance = 1e-10)
})

test_that("params.empirical_dist returns NULL", {
  e <- empirical_dist(1:5)

  expect_null(params(e))
})

test_that("nparams.empirical_dist returns 0", {
  e <- empirical_dist(1:5)

  expect_equal(nparams(e), 0)
})

test_that("sampler.empirical_dist returns a function that resamples from data", {
  data <- c(1, 2, 3, 4, 5)
  e <- empirical_dist(data)
  samp_fn <- sampler(e)

  # The sampler should return a function
  expect_type(samp_fn, "closure")

  # Samples should be from the original data
  set.seed(42)
  samples <- samp_fn(100)
  expect_true(all(samples %in% data))
})

test_that("sampler.empirical_dist generates correct number of samples", {
  e <- empirical_dist(1:10)
  samp_fn <- sampler(e)

  samples <- samp_fn(50)
  expect_length(samples, 50)
})

test_that("sampler.empirical_dist works for multivariate data", {
  data <- matrix(1:10, nrow = 5, ncol = 2)
  e <- empirical_dist(data)
  samp_fn <- sampler(e)

  samples <- samp_fn(100)
  expect_true(is.matrix(samples))
  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 2)
})

test_that("density.empirical_dist returns correct probabilities", {
  data <- c(1, 1, 2, 3, 3, 3)
  e <- empirical_dist(data)
  pdf <- density(e)

  # 1 appears 2 times out of 6
  expect_equal(pdf(1), 2 / 6)
  # 2 appears 1 time out of 6
  expect_equal(pdf(2), 1 / 6)
  # 3 appears 3 times out of 6
  expect_equal(pdf(3), 3 / 6)
  # 4 never appears
  expect_equal(pdf(4), 0)
})

test_that("density.empirical_dist handles log argument correctly", {
  data <- c(1, 1, 2)
  e <- empirical_dist(data)
  pdf <- density(e)

  expect_equal(pdf(1, log = TRUE), log(2 / 3), tolerance = 1e-10)
})

test_that("cdf.empirical_dist returns empirical CDF function", {
  data <- c(1, 2, 3, 4, 5)
  e <- empirical_dist(data)
  cdf_fn <- cdf(e)

  # Check CDF values at data points
  expect_equal(cdf_fn(1), 0.2)
  expect_equal(cdf_fn(3), 0.6)
  expect_equal(cdf_fn(5), 1.0)
  expect_equal(cdf_fn(0), 0)
})

test_that("cdf.empirical_dist errors for multivariate distribution", {
  data <- matrix(1:6, nrow = 3, ncol = 2)
  e <- empirical_dist(data)

  expect_error(cdf(e))
})

test_that("sup.empirical_dist returns finite_set of observed values", {
  data <- c(1, 2, 3)
  e <- empirical_dist(data)
  s <- sup(e)

  expect_s3_class(s, "finite_set")
  expect_true(s$has(1))
  expect_true(s$has(2))
  expect_true(s$has(3))
  expect_false(s$has(4))
})

test_that("marginal.empirical_dist returns correct marginal distribution", {
  data <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, ncol = 2)
  e <- empirical_dist(data)

  # Marginal of first column
  marg <- marginal(e, 1)

  expect_s3_class(marg, "empirical_dist")
  expect_equal(dim(marg), 1)
  expect_equal(as.vector(obs(marg)), c(1, 2, 3))
})

test_that("marginal.empirical_dist validates indices", {
  e <- empirical_dist(matrix(1:6, nrow = 3, ncol = 2))

  expect_error(marginal(e, 0))
  expect_error(marginal(e, 3))
})

test_that("conditional.empirical_dist filters data by predicate", {
  data <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  e <- empirical_dist(data)

  # Condition on values > 5
  cond_e <- conditional(e, function(x) x > 5)

  expect_s3_class(cond_e, "empirical_dist")
  expect_equal(nobs(cond_e), 5)
  expect_true(all(obs(cond_e) > 5))
})

test_that("conditional.empirical_dist works with multivariate data", {
  data <- matrix(c(1, 2, 3, 4, 10, 20, 30, 40), nrow = 4, ncol = 2)
  e <- empirical_dist(data)

  # Condition on first column > 2
  cond_e <- conditional(e, function(x) x[1] > 2)

  expect_equal(nobs(cond_e), 2)
})

test_that("rmap.empirical_dist applies function to observations", {
  data <- c(1, 2, 3, 4, 5)
  e <- empirical_dist(data)

  # Square each observation
  mapped <- rmap(e, function(x) x^2)

  expect_s3_class(mapped, "empirical_dist")
  expect_equal(as.vector(obs(mapped)), c(1, 4, 9, 16, 25))
})

test_that("rmap.empirical_dist handles dimension changes", {
  data <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, ncol = 2)
  e <- empirical_dist(data)

  # Sum across columns (reduces to univariate)
  mapped <- rmap(e, function(x) sum(x))

  expect_s3_class(mapped, "empirical_dist")
  expect_equal(dim(mapped), 1)
  expect_equal(as.vector(obs(mapped)), c(11, 22, 33))
})

test_that("expectation.empirical_dist computes correct expectation", {
  data <- c(1, 2, 3, 4, 5)
  e <- empirical_dist(data)

  # E[X] = mean
  expect_equal(expectation(e), mean(data))

  # E[X^2]
  expect_equal(expectation(e, function(x) x^2), mean(data^2))
})

test_that("expectation.empirical_dist computes stats when requested", {
  data <- 1:100
  e <- empirical_dist(data)

  result <- expectation(e, control = list(compute_stats = TRUE))

  expect_type(result, "list")
  expect_true("value" %in% names(result))
  expect_true("ci" %in% names(result))
})
