test_that("realize.dist returns an empirical_dist", {
  d <- normal(mu = 0, var = 1)
  set.seed(42)
  r <- realize(d, n = 100)
  expect_s3_class(r, "empirical_dist")
  expect_true(is_empirical_dist(r))
})

test_that("realize.empirical_dist returns self (identical object)", {
  ed <- empirical_dist(rnorm(50))
  r <- realize(ed)
  expect_identical(r, ed)
})

test_that("realize with default n=10000 gives correct number of observations", {
  d <- exponential(rate = 2)
  set.seed(1)
  r <- realize(d)
  expect_equal(nobs(r), 10000L)
})

test_that("realize with custom n works", {
  d <- normal(mu = 5, var = 2)
  set.seed(1)
  r <- realize(d, n = 500)
  expect_equal(nobs(r), 500L)
})

test_that("realize with n=0 errors", {
  d <- normal()
  expect_error(realize(d, n = 0), "'n' must be a positive integer")
})

test_that("realize with n=-1 errors", {
  d <- normal()
  expect_error(realize(d, n = -1), "'n' must be a positive integer")
})

test_that("realize works on edist objects", {
  x <- normal(mu = 0, var = 1)
  y <- exponential(rate = 1)
  e <- edist(quote(x + y), list(x = x, y = y))

  set.seed(42)
  r <- realize(e, n = 200)
  expect_s3_class(r, "empirical_dist")
  expect_equal(nobs(r), 200L)
})

test_that("cdf.edist auto-fallback works", {
  x <- normal(mu = 0, var = 1)
  y <- exponential(rate = 1)
  e <- edist(quote(x + y), list(x = x, y = y))

  set.seed(42)
  cdf_fn <- cdf(e)
  expect_true(is.function(cdf_fn))

  vals <- cdf_fn(c(-10, 0, 5, 100))
  expect_true(all(vals >= 0 & vals <= 1))
  # extreme values should be near 0 and 1

  expect_true(vals[1] < 0.01)
  expect_true(vals[4] > 0.99)
})

test_that("density.edist auto-fallback works", {
  x <- normal(mu = 0, var = 1)
  y <- exponential(rate = 1)
  e <- edist(quote(x + y), list(x = x, y = y))

  set.seed(42)
  dens_fn <- density(e)
  expect_true(is.function(dens_fn))

  # density at observed values should be non-negative
  val <- dens_fn(1)
  expect_true(is.numeric(val))
  expect_true(val >= 0)
})

test_that("sup.edist auto-fallback works", {
  x <- normal(mu = 0, var = 1)
  y <- exponential(rate = 1)
  e <- edist(quote(x + y), list(x = x, y = y))

  set.seed(42)
  s <- sup(e)
  # sup of a realized empirical_dist is a finite_set
  expect_true(inherits(s, "finite_set"))
})

test_that("conditional.edist auto-fallback works", {
  x <- normal(mu = 5, var = 1)
  y <- exponential(rate = 1)
  e <- edist(quote(x + y), list(x = x, y = y))

  set.seed(42)
  # condition on values > 5 (most should satisfy this since mean ~= 6)
  cond <- conditional(e, function(d) d[1] > 5)
  expect_s3_class(cond, "empirical_dist")
  # all observations in the conditional should be > 5
  expect_true(all(obs(cond) > 5))
})

test_that("rmap.edist auto-fallback works", {
  x <- normal(mu = 0, var = 1)
  y <- exponential(rate = 1)
  e <- edist(quote(x + y), list(x = x, y = y))

  set.seed(42)
  mapped <- rmap(e, function(d) d^2)
  expect_s3_class(mapped, "empirical_dist")
  # squared values must be non-negative
  expect_true(all(obs(mapped) >= 0))
})

test_that("inv_cdf.edist auto-fallback works", {
  x <- normal(mu = 0, var = 1)
  y <- exponential(rate = 1)
  e <- edist(quote(x + y), list(x = x, y = y))

  set.seed(42)
  q_fn <- inv_cdf(e)
  expect_true(is.function(q_fn))

  quantiles <- q_fn(c(0.1, 0.5, 0.9))
  expect_equal(length(quantiles), 3)
  # quantiles should be monotonically non-decreasing
  expect_true(all(diff(quantiles) >= 0))
})

test_that("realize with non-numeric n errors", {
  d <- normal()
  expect_error(realize(d, n = "abc"), "'n' must be a positive integer")
})

test_that("realize preserves distribution properties approximately", {
  d <- normal(mu = 10, var = 4)
  set.seed(123)
  r <- realize(d, n = 50000)
  # mean should be close to 10

  expect_equal(mean(r), 10, tolerance = 0.1)
})

test_that("inv_cdf.empirical_dist works directly", {
  ed <- empirical_dist(c(1, 2, 3, 4, 5))
  q_fn <- inv_cdf(ed)
  expect_true(is.function(q_fn))
  expect_equal(q_fn(0), 1)
  expect_equal(q_fn(1), 5)
  expect_equal(q_fn(0.5), 3)
})
