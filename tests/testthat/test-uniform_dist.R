# Tests for the uniform distribution class

# -- Construction -----------------------------------------------------------

test_that("uniform_dist constructor creates valid object with default params", {
  u <- uniform_dist()

  expect_s3_class(u, "uniform_dist")
  expect_s3_class(u, "univariate_dist")
  expect_s3_class(u, "continuous_dist")
  expect_s3_class(u, "dist")
  expect_equal(u$min, 0)
  expect_equal(u$max, 1)
})

test_that("uniform_dist constructor creates valid object with custom params", {
  u <- uniform_dist(min = -5, max = 10)

  expect_s3_class(u, "uniform_dist")
  expect_equal(u$min, -5)
  expect_equal(u$max, 10)
})

test_that("uniform_dist errors when min >= max", {
  expect_error(uniform_dist(min = 5, max = 5), "'min' must be less than 'max'")
  expect_error(uniform_dist(min = 10, max = 5), "'min' must be less than 'max'")
})

test_that("uniform_dist errors with non-numeric min", {
  expect_error(uniform_dist(min = "a"), "'min' must be a numeric scalar")
  expect_error(uniform_dist(min = c(1, 2)), "'min' must be a numeric scalar")
})

test_that("uniform_dist errors with non-numeric max", {
  expect_error(uniform_dist(max = "b"), "'max' must be a numeric scalar")
  expect_error(uniform_dist(max = c(1, 2)), "'max' must be a numeric scalar")
})

test_that("uniform_dist errors with NA parameters", {
  expect_error(uniform_dist(min = NA), "'min' must be a numeric scalar")
  expect_error(uniform_dist(max = NA), "'max' must be a numeric scalar")
  expect_error(uniform_dist(min = NA_real_), "'min' must be a numeric scalar")
})

# -- is_uniform_dist -------------------------------------------------------

test_that("is_uniform_dist identifies uniform_dist objects correctly", {
  u <- uniform_dist()

  expect_true(is_uniform_dist(u))
  expect_false(is_uniform_dist(list(min = 0, max = 1)))
  expect_false(is_uniform_dist(normal()))
})

# -- Class hierarchy -------------------------------------------------------

test_that("uniform_dist has correct class hierarchy", {
  u <- uniform_dist(min = 2, max = 8)

  classes <- class(u)
  expect_equal(classes, c("uniform_dist", "univariate_dist",
                          "continuous_dist", "dist"))
})

# -- params -----------------------------------------------------------------

test_that("params.uniform_dist returns named vector of parameters", {
  u <- uniform_dist(min = 3, max = 7)
  p <- params(u)

  expect_named(p, c("min", "max"))
  expect_equal(p["min"], c(min = 3))
  expect_equal(p["max"], c(max = 7))
})

# -- mean -------------------------------------------------------------------

test_that("mean.uniform_dist returns (min + max) / 2", {
  expect_equal(mean(uniform_dist()), 0.5)
  expect_equal(mean(uniform_dist(min = 2, max = 8)), 5)
  expect_equal(mean(uniform_dist(min = -10, max = 10)), 0)
})

# -- vcov -------------------------------------------------------------------

test_that("vcov.uniform_dist returns (max - min)^2 / 12", {
  expect_equal(vcov(uniform_dist()), 1 / 12)
  expect_equal(vcov(uniform_dist(min = 0, max = 6)), 36 / 12)
  expect_equal(vcov(uniform_dist(min = -3, max = 3)), 36 / 12)
})

# -- dim --------------------------------------------------------------------

test_that("dim.uniform_dist returns 1", {
  u <- uniform_dist()
  expect_equal(dim(u), 1)
})

# -- format and print -------------------------------------------------------

test_that("format.uniform_dist returns descriptive string", {
  u <- uniform_dist(min = 2, max = 5)
  fmt <- format(u)

  expect_type(fmt, "character")
  expect_match(fmt, "Uniform distribution")
  expect_match(fmt, "min = 2")
  expect_match(fmt, "max = 5")
})

test_that("print.uniform_dist produces output and returns invisibly", {
  u <- uniform_dist(min = 0, max = 1)

  expect_output(print(u), "Uniform distribution")
  expect_output(print(u), "min = 0")
  expect_output(print(u), "max = 1")

  # Returns invisibly
  expect_invisible(print(u))
})

# -- sampler ----------------------------------------------------------------

test_that("sampler.uniform_dist returns a function that generates samples", {
  u <- uniform_dist(min = 2, max = 5)
  samp_fn <- sampler(u)

  expect_type(samp_fn, "closure")

  samples <- samp_fn(100)
  expect_length(samples, 100)
})

test_that("sampler.uniform_dist produces samples within [min, max]", {
  set.seed(42)
  u <- uniform_dist(min = 2, max = 5)
  samples <- sampler(u)(10000)

  expect_true(all(samples >= 2))
  expect_true(all(samples <= 5))
})

test_that("sampler.uniform_dist produces samples with approximately correct mean", {
  set.seed(42)
  u <- uniform_dist(min = 0, max = 10)
  samples <- sampler(u)(10000)

  sample_mean <- sum(samples) / length(samples)
  expect_equal(sample_mean, 5, tolerance = 0.2)
})

# -- density ----------------------------------------------------------------

test_that("density.uniform_dist returns correct probability density", {
  u <- uniform_dist(min = 2, max = 5)
  pdf <- density(u)

  # Inside the interval: 1/(max-min) = 1/3
  expect_equal(pdf(3), dunif(3, min = 2, max = 5), tolerance = 1e-10)
  expect_equal(pdf(2), dunif(2, min = 2, max = 5), tolerance = 1e-10)
  expect_equal(pdf(5), dunif(5, min = 2, max = 5), tolerance = 1e-10)
})

test_that("density.uniform_dist returns zero outside support", {
  u <- uniform_dist(min = 2, max = 5)
  pdf <- density(u)

  expect_equal(pdf(1), 0)
  expect_equal(pdf(6), 0)
})

test_that("density.uniform_dist handles log argument correctly", {
  u <- uniform_dist(min = 2, max = 5)
  pdf <- density(u)

  expect_equal(pdf(3, log = TRUE),
               dunif(3, min = 2, max = 5, log = TRUE),
               tolerance = 1e-10)
})

# -- cdf --------------------------------------------------------------------

test_that("cdf.uniform_dist returns correct cumulative distribution", {
  u <- uniform_dist(min = 0, max = 10)
  cdf_fn <- cdf(u)

  expect_equal(cdf_fn(0), punif(0, min = 0, max = 10), tolerance = 1e-10)
  expect_equal(cdf_fn(5), punif(5, min = 0, max = 10), tolerance = 1e-10)
  expect_equal(cdf_fn(10), punif(10, min = 0, max = 10), tolerance = 1e-10)
  expect_equal(cdf_fn(-1), 0, tolerance = 1e-10)
  expect_equal(cdf_fn(11), 1, tolerance = 1e-10)
})

test_that("cdf.uniform_dist handles log.p argument correctly", {
  u <- uniform_dist(min = 0, max = 1)
  cdf_fn <- cdf(u)

  expect_equal(cdf_fn(0.5, log.p = TRUE),
               punif(0.5, min = 0, max = 1, log.p = TRUE),
               tolerance = 1e-10)
})

# -- inv_cdf ----------------------------------------------------------------

test_that("inv_cdf.uniform_dist returns correct quantiles", {
  u <- uniform_dist(min = 2, max = 8)
  qf <- inv_cdf(u)

  expect_equal(qf(0), 2, tolerance = 1e-10)
  expect_equal(qf(0.5), 5, tolerance = 1e-10)
  expect_equal(qf(1), 8, tolerance = 1e-10)

  # Compare with qunif
  expect_equal(qf(0.25), qunif(0.25, min = 2, max = 8), tolerance = 1e-10)
})

test_that("inv_cdf round-trips with cdf", {
  u <- uniform_dist(min = -3, max = 7)
  cdf_fn <- cdf(u)
  qf <- inv_cdf(u)

  # cdf(inv_cdf(p)) == p
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (p in probs) {
    expect_equal(cdf_fn(qf(p)), p, tolerance = 1e-10)
  }

  # inv_cdf(cdf(x)) == x for x in (min, max)
  vals <- c(-2, 0, 3, 5)
  for (v in vals) {
    expect_equal(qf(cdf_fn(v)), v, tolerance = 1e-10)
  }
})

# -- sup --------------------------------------------------------------------

test_that("sup.uniform_dist returns closed interval [min, max]", {
  u <- uniform_dist(min = 2, max = 5)
  s <- sup(u)

  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), 2)
  expect_equal(s$supremum(), 5)
  expect_true(s$lower_closed)
  expect_true(s$upper_closed)
})

test_that("sup.uniform_dist support contains interior and boundary points", {
  u <- uniform_dist(min = 0, max = 1)
  s <- sup(u)

  expect_true(s$has(0))
  expect_true(s$has(0.5))
  expect_true(s$has(1))
  expect_false(s$has(-0.1))
  expect_false(s$has(1.1))
})
