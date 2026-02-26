# Tests for the chi-squared distribution class

test_that("chi_squared constructor creates valid object", {
  x <- chi_squared(df = 5)

  expect_s3_class(x, "chi_squared")
  expect_s3_class(x, "univariate_dist")
  expect_s3_class(x, "continuous_dist")
  expect_s3_class(x, "dist")
  expect_equal(x$df, 5)
})

test_that("chi_squared constructor rejects invalid df", {
  expect_error(chi_squared(df = 0), "'df' must be a positive scalar")
  expect_error(chi_squared(df = -1), "'df' must be a positive scalar")
  expect_error(chi_squared(df = "a"), "'df' must be a positive scalar")
  expect_error(chi_squared(df = c(1, 2)), "'df' must be a positive scalar")
  expect_error(chi_squared(df = NA_real_), "'df' must be a positive scalar")
})

test_that("is_chi_squared identifies chi_squared objects correctly", {
  x <- chi_squared(df = 3)

  expect_true(is_chi_squared(x))
  expect_false(is_chi_squared(list(df = 3)))
  expect_false(is_chi_squared(normal()))
  expect_false(is_chi_squared(exponential(rate = 1)))
})

test_that("params.chi_squared returns named vector of parameters", {
  x <- chi_squared(df = 7)
  p <- params(x)

  expect_named(p, "df")
  expect_equal(p["df"], c(df = 7))
})

test_that("mean.chi_squared returns df", {
  expect_equal(mean(chi_squared(df = 3)), 3)
  expect_equal(mean(chi_squared(df = 10)), 10)
  expect_equal(mean(chi_squared(df = 0.5)), 0.5)
})

test_that("vcov.chi_squared returns 2*df", {
  expect_equal(vcov(chi_squared(df = 3)), 6)
  expect_equal(vcov(chi_squared(df = 10)), 20)
  expect_equal(vcov(chi_squared(df = 0.5)), 1)
})

test_that("dim.chi_squared returns 1 for univariate distribution", {
  expect_equal(dim(chi_squared(df = 5)), 1)
})

test_that("format.chi_squared returns descriptive string", {
  x <- chi_squared(df = 4)
  expect_equal(format(x), "Chi-squared distribution (df = 4)")
})

test_that("print.chi_squared produces output and returns invisibly", {
  x <- chi_squared(df = 4)

  expect_output(print(x), "Chi-squared distribution")
  expect_output(print(x), "df = 4")
  expect_invisible(print(x))
})

test_that("sampler.chi_squared returns a function that generates n samples", {
  x <- chi_squared(df = 5)
  samp_fn <- sampler(x)

  expect_type(samp_fn, "closure")

  samples <- samp_fn(100)
  expect_length(samples, 100)

  # Samples should be non-negative (chi-squared is on (0, Inf))
  expect_true(all(samples >= 0))
})

test_that("sampler.chi_squared produces samples with approximately correct mean", {
  set.seed(42)
  x <- chi_squared(df = 5)
  samp_fn <- sampler(x)
  samples <- samp_fn(10000)

  # Mean should be df = 5
  sample_mean <- sum(samples) / length(samples)
  expect_equal(sample_mean, 5, tolerance = 0.2)
})

test_that("density.chi_squared matches dchisq", {
  x <- chi_squared(df = 3)
  pdf <- density(x)

  # Compare with dchisq at known points
  expect_equal(pdf(1), dchisq(1, df = 3), tolerance = 1e-10)
  expect_equal(pdf(2), dchisq(2, df = 3), tolerance = 1e-10)
  expect_equal(pdf(5), dchisq(5, df = 3), tolerance = 1e-10)
})

test_that("density.chi_squared handles log argument correctly", {
  x <- chi_squared(df = 4)
  pdf <- density(x)

  expect_equal(pdf(2, log = TRUE), dchisq(2, df = 4, log = TRUE),
               tolerance = 1e-10)
})

test_that("density.chi_squared returns zero for negative values", {
  x <- chi_squared(df = 3)
  pdf <- density(x)

  expect_equal(pdf(-1), 0)
})

test_that("cdf.chi_squared matches pchisq", {
  x <- chi_squared(df = 4)
  cdf_fn <- cdf(x)

  expect_equal(cdf_fn(0), pchisq(0, df = 4), tolerance = 1e-10)
  expect_equal(cdf_fn(1), pchisq(1, df = 4), tolerance = 1e-10)
  expect_equal(cdf_fn(5), pchisq(5, df = 4), tolerance = 1e-10)
})

test_that("cdf.chi_squared handles log.p argument correctly", {
  x <- chi_squared(df = 3)
  cdf_fn <- cdf(x)

  expect_equal(cdf_fn(2, log.p = TRUE), pchisq(2, df = 3, log.p = TRUE),
               tolerance = 1e-10)
})

test_that("inv_cdf.chi_squared round-trips with cdf", {
  x <- chi_squared(df = 5)
  cdf_fn <- cdf(x)
  qf <- inv_cdf(x)

  # Round-trip: inv_cdf(cdf(t)) == t
  for (t in c(0.5, 1, 2, 5, 10)) {
    p <- cdf_fn(t)
    expect_equal(qf(p), t, tolerance = 1e-10)
  }

  # Round-trip: cdf(inv_cdf(p)) == p
  for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    t <- qf(p)
    expect_equal(cdf_fn(t), p, tolerance = 1e-10)
  }
})

test_that("inv_cdf.chi_squared matches qchisq", {
  x <- chi_squared(df = 6)
  qf <- inv_cdf(x)

  expect_equal(qf(0.5), qchisq(0.5, df = 6), tolerance = 1e-10)
  expect_equal(qf(0.95), qchisq(0.95, df = 6), tolerance = 1e-10)
})

test_that("sup.chi_squared returns interval (0, Inf)", {
  x <- chi_squared(df = 3)
  s <- sup(x)

  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), 0)
  expect_equal(s$supremum(), Inf)
  expect_false(s$lower_closed)
})

test_that("surv.chi_squared equals 1 - cdf at test points", {
  x <- chi_squared(df = 4)
  cdf_fn <- cdf(x)
  sf <- surv(x)

  for (t in c(0.5, 1, 2, 5, 10)) {
    expect_equal(sf(t), 1 - cdf_fn(t), tolerance = 1e-10)
  }
})

test_that("surv.chi_squared matches pchisq lower.tail=FALSE", {
  x <- chi_squared(df = 3)
  sf <- surv(x)

  expect_equal(sf(2), pchisq(2, df = 3, lower.tail = FALSE), tolerance = 1e-10)
  expect_equal(sf(5), pchisq(5, df = 3, lower.tail = FALSE), tolerance = 1e-10)
})

test_that("hazard.chi_squared equals f(t)/S(t) for several test points", {
  x <- chi_squared(df = 4)
  h <- hazard(x)
  pdf <- density(x)
  sf <- surv(x)

  for (t in c(1, 2, 5)) {
    expected_h <- pdf(t) / sf(t)
    expect_equal(h(t), expected_h, tolerance = 1e-10)
  }
})

test_that("hazard.chi_squared handles log argument correctly", {
  x <- chi_squared(df = 4)
  h <- hazard(x)
  pdf <- density(x)
  sf <- surv(x)

  expected_h <- pdf(2) / sf(2)
  expect_equal(h(2, log = TRUE), log(expected_h), tolerance = 1e-10)
})

test_that("chi_squared works with non-integer df", {
  x <- chi_squared(df = 2.5)

  expect_equal(mean(x), 2.5)
  expect_equal(vcov(x), 5)

  pdf <- density(x)
  expect_equal(pdf(1), dchisq(1, df = 2.5), tolerance = 1e-10)

  cdf_fn <- cdf(x)
  expect_equal(cdf_fn(1), pchisq(1, df = 2.5), tolerance = 1e-10)
})
