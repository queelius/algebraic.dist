# Tests for the exponential distribution class

test_that("exponential constructor creates valid object", {
  # Given: A rate parameter
  rate <- 2

  # When: Creating an exponential distribution
  e <- exponential(rate = rate)

  # Then: Object has correct class hierarchy and parameter
  expect_s3_class(e, "exponential")
  expect_s3_class(e, "univariate_dist")
  expect_s3_class(e, "continuous_dist")
  expect_s3_class(e, "dist")
  expect_equal(e$rate, 2)
})

test_that("is_exponential identifies exponential objects correctly", {
  e <- exponential(rate = 1)

  expect_true(is_exponential(e))
  expect_false(is_exponential(list(rate = 1)))
  expect_false(is_exponential(normal()))
})

test_that("mean.exponential returns correct mean", {
  # Given: An exponential distribution with rate = 2
  e <- exponential(rate = 2)

  # When/Then: Mean should be 1/rate = 0.5
  expect_equal(mean(e), 0.5)
})

test_that("params.exponential returns named vector of parameters", {
  e <- exponential(rate = 3)
  p <- params(e)

  expect_named(p, "rate")
  expect_equal(p["rate"], c(rate = 3))
})

test_that("dim.exponential returns 1 for univariate distribution", {
  e <- exponential(rate = 1)

  expect_equal(dim(e), 1)
})

test_that("sampler.exponential returns a function that generates samples", {
  e <- exponential(rate = 1)
  samp_fn <- sampler(e)

  # The sampler should return a function
  expect_type(samp_fn, "closure")

  # The function should generate the requested number of samples
  samples <- samp_fn(100)
  expect_length(samples, 100)

  # Samples should be non-negative (exponential is on (0, Inf))
  expect_true(all(samples >= 0))
})

test_that("sampler.exponential produces samples with approximately correct mean", {
  set.seed(42)
  e <- exponential(rate = 2)
  samp_fn <- sampler(e)
  samples <- samp_fn(10000)

  # Mean should be 1/rate = 0.5
  # Use sum/length to compute mean (avoid mean.default override in package)
  sample_mean <- sum(samples) / length(samples)
  expect_type(sample_mean, "double")
  expect_equal(sample_mean, 0.5, tolerance = 0.1)
})

test_that("density.exponential returns correct probability density", {
  e <- exponential(rate = 2)
  pdf <- density(e)

  # Compare with dexp at known points
  expect_equal(pdf(0), dexp(0, rate = 2), tolerance = 1e-10)
  expect_equal(pdf(1), dexp(1, rate = 2), tolerance = 1e-10)
  expect_equal(pdf(2), dexp(2, rate = 2), tolerance = 1e-10)
})

test_that("density.exponential handles log argument correctly", {
  e <- exponential(rate = 2)
  pdf <- density(e)

  expect_equal(pdf(1, log = TRUE), dexp(1, rate = 2, log = TRUE), tolerance = 1e-10)
})

test_that("density.exponential returns zero for negative values", {
  e <- exponential(rate = 1)
  pdf <- density(e)

  expect_equal(pdf(-1), 0)
})

test_that("cdf.exponential returns correct cumulative distribution", {
  e <- exponential(rate = 2)
  cdf_fn <- cdf(e)

  # Compare with pexp at known points
  expect_equal(cdf_fn(0), pexp(0, rate = 2), tolerance = 1e-10)
  expect_equal(cdf_fn(1), pexp(1, rate = 2), tolerance = 1e-10)
  expect_equal(cdf_fn(0.5), pexp(0.5, rate = 2), tolerance = 1e-10)
})

test_that("cdf.exponential handles log.p argument correctly", {
  e <- exponential(rate = 1)
  cdf_fn <- cdf(e)

  expect_equal(cdf_fn(1, log.p = TRUE), pexp(1, rate = 1, log.p = TRUE), tolerance = 1e-10)
})

test_that("inv_cdf.exponential returns correct quantiles", {
  e <- exponential(rate = 2)
  qf <- inv_cdf(e)

  # Compare with qexp at known points
  expect_equal(qf(0.5), qexp(0.5, rate = 2), tolerance = 1e-10)
  expect_equal(qf(0.95), qexp(0.95, rate = 2), tolerance = 1e-10)
})

test_that("surv.exponential returns correct survival function", {
  e <- exponential(rate = 2)
  sf <- surv(e)

  # Survival function S(t) = 1 - F(t) = P(X > t)
  expect_equal(sf(0), 1, tolerance = 1e-10)
  expect_equal(sf(1), pexp(1, rate = 2, lower.tail = FALSE), tolerance = 1e-10)
})

test_that("hazard.exponential returns constant hazard rate", {
  e <- exponential(rate = 2)
  h <- hazard(e)

  # For exponential distribution, hazard function is constant = rate
  expect_equal(h(0.5), 2)
  expect_equal(h(10), 2)
})

test_that("hazard.exponential returns zero for negative values", {
  e <- exponential(rate = 2)
  h <- hazard(e)

  expect_equal(h(-1), 0)
})

test_that("hazard.exponential handles log argument correctly", {
  e <- exponential(rate = 2)
  h <- hazard(e)

  expect_equal(h(1, log = TRUE), log(2))
  expect_equal(h(-1, log = TRUE), -Inf)
})

test_that("sup.exponential returns correct support interval", {
  e <- exponential(rate = 1)
  s <- sup(e)

  # Exponential distribution has support (0, Inf)
  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), 0)
  expect_equal(s$supremum(), Inf)
  # Lower bound should be open (not including 0)
  expect_false(s$lower_closed)
})

test_that("print.exponential produces output without error", {
  e <- exponential(rate = 2)

  expect_output(print(e), "Exponential distribution")
  expect_output(print(e), "rate = 2")
})
