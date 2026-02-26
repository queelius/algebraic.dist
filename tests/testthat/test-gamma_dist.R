# Tests for the gamma distribution class

# --- Construction ---

test_that("gamma_dist constructor creates valid object with correct params", {
  g <- gamma_dist(shape = 2, rate = 3)

  expect_s3_class(g, "gamma_dist")
  expect_s3_class(g, "univariate_dist")
  expect_s3_class(g, "continuous_dist")
  expect_s3_class(g, "dist")
  expect_equal(g$shape, 2)
  expect_equal(g$rate, 3)
})

test_that("gamma_dist constructor rejects invalid shape", {
  expect_error(gamma_dist(shape = -1, rate = 1), "'shape' must be a positive scalar")
  expect_error(gamma_dist(shape = 0, rate = 1), "'shape' must be a positive scalar")
  expect_error(gamma_dist(shape = "a", rate = 1), "'shape' must be a positive scalar")
  expect_error(gamma_dist(shape = c(1, 2), rate = 1), "'shape' must be a positive scalar")
  expect_error(gamma_dist(shape = NA_real_, rate = 1), "'shape' must be a positive scalar")
})

test_that("gamma_dist constructor rejects invalid rate", {
  expect_error(gamma_dist(shape = 1, rate = -1), "'rate' must be a positive scalar")
  expect_error(gamma_dist(shape = 1, rate = 0), "'rate' must be a positive scalar")
  expect_error(gamma_dist(shape = 1, rate = "b"), "'rate' must be a positive scalar")
  expect_error(gamma_dist(shape = 1, rate = c(1, 2)), "'rate' must be a positive scalar")
  expect_error(gamma_dist(shape = 1, rate = NA_real_), "'rate' must be a positive scalar")
})

# --- is_gamma_dist ---

test_that("is_gamma_dist identifies gamma_dist objects correctly", {
  g <- gamma_dist(shape = 2, rate = 1)

  expect_true(is_gamma_dist(g))
  expect_false(is_gamma_dist(list(shape = 2, rate = 1)))
  expect_false(is_gamma_dist(normal()))
  expect_false(is_gamma_dist(exponential(rate = 1)))
})

# --- params ---

test_that("params.gamma_dist returns named vector of parameters", {
  g <- gamma_dist(shape = 2.5, rate = 0.7)
  p <- params(g)

  expect_named(p, c("shape", "rate"))
  expect_equal(p["shape"], c(shape = 2.5))
  expect_equal(p["rate"], c(rate = 0.7))
})

# --- mean ---

test_that("mean.gamma_dist returns shape/rate", {
  g <- gamma_dist(shape = 3, rate = 2)
  expect_equal(mean(g), 1.5)

  g2 <- gamma_dist(shape = 5, rate = 10)
  expect_equal(mean(g2), 0.5)
})

# --- vcov ---

test_that("vcov.gamma_dist returns shape/rate^2", {
  g <- gamma_dist(shape = 3, rate = 2)
  expect_equal(vcov(g), 0.75)

  g2 <- gamma_dist(shape = 4, rate = 1)
  expect_equal(vcov(g2), 4)
})

# --- dim ---

test_that("dim.gamma_dist returns 1 for univariate distribution", {
  g <- gamma_dist(shape = 1, rate = 1)
  expect_equal(dim(g), 1)
})

# --- format / print ---

test_that("format.gamma_dist returns correct string", {
  g <- gamma_dist(shape = 2, rate = 3)
  expect_equal(format(g), "Gamma distribution (shape = 2, rate = 3)")
})

test_that("print.gamma_dist outputs to console", {
  g <- gamma_dist(shape = 2, rate = 3)
  expect_output(print(g), "Gamma distribution")
  expect_output(print(g), "shape = 2")
  expect_output(print(g), "rate = 3")
})

test_that("print.gamma_dist returns object invisibly", {
  g <- gamma_dist(shape = 1, rate = 1)
  out <- capture.output(ret <- print(g))
  expect_identical(ret, g)
})

# --- sampler ---

test_that("sampler.gamma_dist returns a function that generates samples", {
  g <- gamma_dist(shape = 2, rate = 1)
  samp_fn <- sampler(g)

  expect_type(samp_fn, "closure")

  samples <- samp_fn(100)
  expect_length(samples, 100)

  # Samples should be non-negative (gamma is on (0, Inf))
  expect_true(all(samples >= 0))
})

test_that("sampler.gamma_dist produces samples with approximately correct mean", {
  set.seed(42)
  g <- gamma_dist(shape = 3, rate = 2)
  samp_fn <- sampler(g)
  samples <- samp_fn(10000)

  # Mean should be shape/rate = 1.5
  sample_mean <- sum(samples) / length(samples)
  expect_equal(sample_mean, 1.5, tolerance = 0.1)
})

# --- density ---

test_that("density.gamma_dist returns correct probability density", {
  g <- gamma_dist(shape = 2, rate = 3)
  pdf <- density(g)

  # Compare with dgamma at known points
  expect_equal(pdf(0.5), dgamma(0.5, shape = 2, rate = 3), tolerance = 1e-10)
  expect_equal(pdf(1), dgamma(1, shape = 2, rate = 3), tolerance = 1e-10)
  expect_equal(pdf(2), dgamma(2, shape = 2, rate = 3), tolerance = 1e-10)
})

test_that("density.gamma_dist handles log argument correctly", {
  g <- gamma_dist(shape = 2, rate = 3)
  pdf <- density(g)

  expect_equal(pdf(1, log = TRUE), dgamma(1, shape = 2, rate = 3, log = TRUE),
               tolerance = 1e-10)
})

test_that("density.gamma_dist returns zero for negative values", {
  g <- gamma_dist(shape = 2, rate = 1)
  pdf <- density(g)

  expect_equal(pdf(-1), 0)
})

# --- cdf ---

test_that("cdf.gamma_dist returns correct cumulative distribution", {
  g <- gamma_dist(shape = 2, rate = 3)
  cdf_fn <- cdf(g)

  expect_equal(cdf_fn(0), pgamma(0, shape = 2, rate = 3), tolerance = 1e-10)
  expect_equal(cdf_fn(1), pgamma(1, shape = 2, rate = 3), tolerance = 1e-10)
  expect_equal(cdf_fn(0.5), pgamma(0.5, shape = 2, rate = 3), tolerance = 1e-10)
})

test_that("cdf.gamma_dist handles log.p argument correctly", {
  g <- gamma_dist(shape = 2, rate = 3)
  cdf_fn <- cdf(g)

  expect_equal(cdf_fn(1, log.p = TRUE),
               pgamma(1, shape = 2, rate = 3, log.p = TRUE),
               tolerance = 1e-10)
})

# --- inv_cdf ---

test_that("inv_cdf.gamma_dist returns correct quantiles", {
  g <- gamma_dist(shape = 2, rate = 3)
  qf <- inv_cdf(g)

  expect_equal(qf(0.5), qgamma(0.5, shape = 2, rate = 3), tolerance = 1e-10)
  expect_equal(qf(0.95), qgamma(0.95, shape = 2, rate = 3), tolerance = 1e-10)
})

test_that("inv_cdf.gamma_dist round-trips with cdf", {
  g <- gamma_dist(shape = 2, rate = 3)
  cdf_fn <- cdf(g)
  qf <- inv_cdf(g)

  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (p in probs) {
    expect_equal(cdf_fn(qf(p)), p, tolerance = 1e-10)
  }
})

# --- sup ---

test_that("sup.gamma_dist returns correct support interval", {
  g <- gamma_dist(shape = 2, rate = 1)
  s <- sup(g)

  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), 0)
  expect_equal(s$supremum(), Inf)
  # Lower bound should be open (not including 0)
  expect_false(s$lower_closed)
})

# --- hazard ---

test_that("hazard.gamma_dist equals f(t)/S(t) for several test points", {
  g <- gamma_dist(shape = 2, rate = 3)
  h <- hazard(g)

  test_points <- c(0.1, 0.5, 1, 2, 5)
  for (t in test_points) {
    expected <- dgamma(t, shape = 2, rate = 3) /
                pgamma(t, shape = 2, rate = 3, lower.tail = FALSE)
    expect_equal(h(t), expected, tolerance = 1e-10,
                 label = paste("hazard at t =", t))
  }
})

test_that("hazard.gamma_dist handles log argument correctly", {
  g <- gamma_dist(shape = 2, rate = 3)
  h <- hazard(g)

  t <- 1
  expected_h <- dgamma(t, shape = 2, rate = 3) /
                pgamma(t, shape = 2, rate = 3, lower.tail = FALSE)
  expect_equal(h(t, log = TRUE), log(expected_h), tolerance = 1e-10)
})

# --- surv ---

test_that("surv.gamma_dist returns 1 - cdf at test points", {
  g <- gamma_dist(shape = 2, rate = 3)
  sf <- surv(g)
  cdf_fn <- cdf(g)

  test_points <- c(0, 0.1, 0.5, 1, 2, 5)
  for (t in test_points) {
    expect_equal(sf(t), 1 - cdf_fn(t), tolerance = 1e-10,
                 label = paste("survival at t =", t))
  }
})

test_that("surv.gamma_dist handles log.p argument", {
  g <- gamma_dist(shape = 2, rate = 3)
  sf <- surv(g)

  expect_equal(sf(1, log.p = TRUE),
               pgamma(1, shape = 2, rate = 3, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-10)
})

test_that("surv.gamma_dist at t=0 is 1", {
  g <- gamma_dist(shape = 2, rate = 3)
  sf <- surv(g)
  expect_equal(sf(0), 1, tolerance = 1e-10)
})

# --- Special case: exponential is gamma(1, rate) ---

test_that("gamma_dist(1, rate) matches exponential distribution", {
  set.seed(123)
  rate <- 2
  g <- gamma_dist(shape = 1, rate = rate)

  # Mean
  expect_equal(mean(g), 1 / rate)

  # Density at several points
  pdf <- density(g)
  for (t in c(0.5, 1, 2)) {
    expect_equal(pdf(t), dexp(t, rate = rate), tolerance = 1e-10)
  }

  # CDF
  cdf_fn <- cdf(g)
  for (t in c(0.5, 1, 2)) {
    expect_equal(cdf_fn(t), pexp(t, rate = rate), tolerance = 1e-10)
  }
})
