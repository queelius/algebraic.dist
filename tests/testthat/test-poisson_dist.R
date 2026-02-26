# Tests for the poisson_dist distribution class

test_that("poisson_dist constructor creates valid object", {
  p <- poisson_dist(lambda = 3)

  expect_s3_class(p, "poisson_dist")
  expect_equal(p$lambda, 3)
})

test_that("poisson_dist constructor rejects invalid lambda", {
  # zero

expect_error(poisson_dist(lambda = 0))
  # negative
  expect_error(poisson_dist(lambda = -1))
  # non-numeric
  expect_error(poisson_dist(lambda = "abc"))
  # vector
  expect_error(poisson_dist(lambda = c(1, 2)))
  # NULL
  expect_error(poisson_dist(lambda = NULL))
  # NA
  expect_error(poisson_dist(lambda = NA), "'lambda' must be a positive scalar")
  expect_error(poisson_dist(lambda = NA_real_), "'lambda' must be a positive scalar")
})

test_that("poisson_dist has correct class hierarchy", {
  p <- poisson_dist(lambda = 5)

  expect_s3_class(p, "poisson_dist")
  expect_s3_class(p, "univariate_dist")
  expect_s3_class(p, "discrete_dist")
  expect_s3_class(p, "dist")
})

test_that("is_poisson_dist identifies poisson objects correctly", {
  p <- poisson_dist(lambda = 2)

  expect_true(is_poisson_dist(p))
  expect_false(is_poisson_dist(list(lambda = 2)))
  expect_false(is_poisson_dist(normal()))
  expect_false(is_poisson_dist(exponential(rate = 1)))
})

test_that("params.poisson_dist returns named vector", {
  p <- poisson_dist(lambda = 4.5)
  par <- params(p)

  expect_named(par, "lambda")
  expect_equal(par[["lambda"]], 4.5)
})

test_that("mean.poisson_dist returns lambda", {
  p <- poisson_dist(lambda = 7)

  expect_equal(mean(p), 7)
})

test_that("vcov.poisson_dist returns lambda", {
  p <- poisson_dist(lambda = 3.5)

  expect_equal(vcov(p), 3.5)
})

test_that("dim.poisson_dist returns 1", {
  p <- poisson_dist(lambda = 1)

  expect_equal(dim(p), 1)
})

test_that("format.poisson_dist produces correct string", {
  p <- poisson_dist(lambda = 2.5)

  expect_match(format(p), "Poisson distribution")
  expect_match(format(p), "lambda = 2.5")
})

test_that("print.poisson_dist produces output without error", {
  p <- poisson_dist(lambda = 2)

  expect_output(print(p), "Poisson distribution")
  expect_output(print(p), "lambda = 2")
})

test_that("print.poisson_dist returns invisibly", {
  p <- poisson_dist(lambda = 1)

  ret <- withVisible(print(p))
  expect_false(ret$visible)
  expect_identical(ret$value, p)
})

test_that("sampler.poisson_dist returns a function that generates integer samples", {
  set.seed(42)
  p <- poisson_dist(lambda = 5)
  samp_fn <- sampler(p)

  expect_type(samp_fn, "closure")

  samples <- samp_fn(1000)
  expect_length(samples, 1000)

  # Samples should be non-negative integers
  expect_true(all(samples >= 0))
  expect_true(all(samples == floor(samples)))
})

test_that("sampler.poisson_dist produces samples with approximately correct mean", {
  set.seed(42)
  p <- poisson_dist(lambda = 10)
  samples <- sampler(p)(10000)

  sample_mean <- sum(samples) / length(samples)
  expect_equal(sample_mean, 10, tolerance = 0.5)
})

test_that("density.poisson_dist matches dpois", {
  p <- poisson_dist(lambda = 3)
  pmf <- density(p)

  for (k in 0:10) {
    expect_equal(pmf(k), dpois(k, lambda = 3), tolerance = 1e-12)
  }
})

test_that("density.poisson_dist handles log argument", {
  p <- poisson_dist(lambda = 4)
  pmf <- density(p)

  for (k in 0:8) {
    expect_equal(pmf(k, log = TRUE), dpois(k, lambda = 4, log = TRUE),
                 tolerance = 1e-12)
  }
})

test_that("density.poisson_dist handles vector input", {
  p <- poisson_dist(lambda = 2)
  pmf <- density(p)

  k <- 0:5
  expect_equal(pmf(k), dpois(k, lambda = 2), tolerance = 1e-12)
})

test_that("cdf.poisson_dist matches ppois", {
  p <- poisson_dist(lambda = 3)
  cdf_fn <- cdf(p)

  for (q in c(0, 1, 2, 5, 10)) {
    expect_equal(cdf_fn(q), ppois(q, lambda = 3), tolerance = 1e-12)
  }
})

test_that("cdf.poisson_dist handles log.p argument", {
  p <- poisson_dist(lambda = 5)
  cdf_fn <- cdf(p)

  expect_equal(cdf_fn(3, log.p = TRUE), ppois(3, lambda = 5, log.p = TRUE),
               tolerance = 1e-12)
})

test_that("inv_cdf.poisson_dist works correctly", {
  p <- poisson_dist(lambda = 4)
  qf <- inv_cdf(p)

  for (prob in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    expect_equal(qf(prob), qpois(prob, lambda = 4))
  }
})

test_that("inv_cdf.poisson_dist handles lower.tail and log.p", {
  p <- poisson_dist(lambda = 3)
  qf <- inv_cdf(p)

  # lower.tail = FALSE
  expect_equal(qf(0.1, lower.tail = FALSE),
               qpois(0.1, lambda = 3, lower.tail = FALSE))

  # log.p = TRUE
  expect_equal(qf(log(0.5), log.p = TRUE),
               qpois(log(0.5), lambda = 3, log.p = TRUE))
})

test_that("sup.poisson_dist returns countable_set with lower bound 0", {
  p <- poisson_dist(lambda = 1)
  s <- sup(p)

  expect_s3_class(s, "countable_set")
  expect_equal(s$lower_bound, 0L)
  expect_equal(infimum(s), 0L)
  expect_equal(supremum(s), Inf)
})

test_that("expectation.poisson_dist with identity gives lambda", {
  p <- poisson_dist(lambda = 7.3)

  result <- expectation(p, function(k) k)
  expect_equal(result, 7.3, tolerance = 1e-10)
})

test_that("expectation.poisson_dist with k^2 gives E[X^2] = lambda + lambda^2", {
  lambda <- 4.2
  p <- poisson_dist(lambda = lambda)

  result <- expectation(p, function(k) k^2)
  expected <- lambda + lambda^2
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("expectation.poisson_dist with constant function gives constant", {
  p <- poisson_dist(lambda = 3)

  result <- expectation(p, function(k) 5)
  expect_equal(result, 5, tolerance = 1e-10)
})

test_that("expectation.poisson_dist variance via E[X^2] - E[X]^2 equals lambda", {
  lambda <- 6.1
  p <- poisson_dist(lambda = lambda)

  ex2 <- expectation(p, function(k) k^2)
  ex <- expectation(p, function(k) k)
  var_result <- ex2 - ex^2
  expect_equal(var_result, lambda, tolerance = 1e-10)
})

test_that("expectation.poisson_dist works for small lambda", {
  p <- poisson_dist(lambda = 0.01)

  result <- expectation(p, function(k) k)
  expect_equal(result, 0.01, tolerance = 1e-9)
})

test_that("expectation.poisson_dist works for large lambda", {
  p <- poisson_dist(lambda = 100)

  result <- expectation(p, function(k) k)
  expect_equal(result, 100, tolerance = 1e-8)
})
