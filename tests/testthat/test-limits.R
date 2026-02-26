# Tests for limiting distribution builder functions (R/limits.R)

# =============================================================================
# clt() — Central Limit Theorem
# =============================================================================

test_that("clt returns normal(0, var) for univariate normal input", {
  # Given: A normal distribution with mu = 5, var = 3
  n <- normal(mu = 5, var = 3)

  # When: Computing the CLT limit
  result <- clt(n)

  # Then: The result is a normal with mu = 0 and var = 3
  expect_true(is_normal(result))
  expect_equal(mean(result), 0)
  expect_equal(vcov(result), 3)
})

test_that("clt returns normal(0, 1/rate^2) for exponential input", {
  # Given: An exponential distribution with rate = 2 (var = 1/4 = 0.25)
  e <- exponential(rate = 2)

  # When: Computing the CLT limit
  result <- clt(e)

  # Then: The result is a normal with mu = 0 and var = 0.25
  expect_true(is_normal(result))
  expect_equal(mean(result), 0)
  expect_equal(vcov(result), 0.25)
})

test_that("clt returns mvn(0, sigma) for multivariate normal input", {
  # Given: A 2D MVN with mean c(1,2) and sigma = I_2
  m <- mvn(mu = c(1, 2), sigma = diag(2))

  # When: Computing the CLT limit
  result <- clt(m)

  # Then: The result is an mvn with mu = c(0,0) and sigma = I_2
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(0, 0))
  expect_equal(vcov(result), diag(2))
})

test_that("clt preserves non-identity covariance matrix", {
  # Given: A 2D MVN with a non-trivial covariance matrix
  sig <- matrix(c(4, 1, 1, 2), 2, 2)
  m <- mvn(mu = c(10, 20), sigma = sig)

  # When: Computing the CLT limit
  result <- clt(m)

  # Then: The result has zero mean and the same covariance
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(0, 0))
  expect_equal(vcov(result), sig)
})

test_that("clt errors on non-dist input", {
  expect_error(clt(42), "'base_dist' must be a 'dist' object")
  expect_error(clt("hello"), "'base_dist' must be a 'dist' object")
  expect_error(clt(list(mu = 0, var = 1)), "'base_dist' must be a 'dist' object")
})

test_that("clt works with gamma distribution", {
  # Given: Gamma(shape=4, rate=2), var = shape/rate^2 = 4/4 = 1

  g <- gamma_dist(shape = 4, rate = 2)

  # When: Computing the CLT limit
  result <- clt(g)

  # Then: normal(0, 1)
  expect_true(is_normal(result))
  expect_equal(mean(result), 0)
  expect_equal(vcov(result), 1)
})

# =============================================================================
# lln() — Law of Large Numbers
# =============================================================================

test_that("lln returns normal(mu, 0) for univariate normal input", {
  # Given: A normal distribution with mu = 5, var = 3
  n <- normal(mu = 5, var = 3)

  # When: Computing the LLN limit
  result <- lln(n)

  # Then: The result is a degenerate normal at mu = 5 with var = 0
  expect_true(is_normal(result))
  expect_equal(mean(result), 5)
  expect_equal(vcov(result), 0)
})

test_that("lln returns normal(1/rate, 0) for exponential input", {
  # Given: An exponential distribution with rate = 2 (mean = 0.5)
  e <- exponential(rate = 2)

  # When: Computing the LLN limit
  result <- lln(e)

  # Then: The result is a degenerate normal at mu = 0.5 with var = 0
  expect_true(is_normal(result))
  expect_equal(mean(result), 0.5)
  expect_equal(vcov(result), 0)
})

test_that("lln returns mvn with zero covariance for multivariate input", {
  # Given: A 2D MVN with mean c(3, 7)
  m <- mvn(mu = c(3, 7), sigma = matrix(c(2, 0.5, 0.5, 3), 2, 2))

  # When: Computing the LLN limit
  result <- lln(m)

  # Then: The result is an mvn with mu = c(3, 7) and zero covariance
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(3, 7))
  expect_equal(vcov(result), matrix(0, 2, 2))
})

test_that("lln errors on non-dist input", {
  expect_error(lln(42), "'base_dist' must be a 'dist' object")
  expect_error(lln("hello"), "'base_dist' must be a 'dist' object")
  expect_error(lln(data.frame(x = 1:5)), "'base_dist' must be a 'dist' object")
})

test_that("lln works with gamma distribution", {
  # Given: Gamma(shape=3, rate=2), mean = 3/2 = 1.5
  g <- gamma_dist(shape = 3, rate = 2)

  # When: Computing the LLN limit
  result <- lln(g)

  # Then: normal(1.5, 0)
  expect_true(is_normal(result))
  expect_equal(mean(result), 1.5)
  expect_equal(vcov(result), 0)
})

# =============================================================================
# delta_clt() — Delta Method CLT
# =============================================================================

test_that("delta_clt returns correct limit for univariate g(x) = x^2", {
  # Given: Exponential(rate=2), mean = 0.5, var = 0.25
  #   g(x) = x^2, dg(x) = 2x
  #   Limiting variance = dg(mu)^2 * var = (2 * 0.5)^2 * 0.25 = 1 * 0.25 = 0.25
  e <- exponential(rate = 2)
  g <- function(x) x^2
  dg <- function(x) 2 * x

  # When: Computing the delta method CLT limit
  result <- delta_clt(e, g, dg)

  # Then: normal(0, 0.25)
  expect_true(is_normal(result))
  expect_equal(mean(result), 0)
  expect_equal(vcov(result), 0.25)
})

test_that("delta_clt returns correct limit for normal input with g(x) = exp(x)", {
  # Given: Normal(mu=1, var=2)
  #   g(x) = exp(x), dg(x) = exp(x)
  #   Limiting variance = exp(1)^2 * 2 = e^2 * 2
  n <- normal(mu = 1, var = 2)
  g <- function(x) exp(x)
  dg <- function(x) exp(x)

  # When: Computing the delta method CLT limit
  result <- delta_clt(n, g, dg)

  # Then: normal(0, exp(1)^2 * 2)
  expect_true(is_normal(result))
  expect_equal(mean(result), 0)
  expect_equal(vcov(result), exp(1)^2 * 2, tolerance = 1e-10)
})

test_that("delta_clt returns correct limit for multivariate input", {
  # Given: MVN(mu = c(1,2), sigma = I_2)
  #   g: R^2 -> R^2 is identity, so J = I_2
  #   Limiting covariance = I_2 %*% I_2 %*% t(I_2) = I_2
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  g <- function(x) x
  dg <- function(x) diag(2)

  # When: Computing the delta method CLT limit
  result <- delta_clt(m, g, dg)

  # Then: mvn(c(0,0), I_2)
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(0, 0))
  expect_equal(vcov(result), diag(2))
})

test_that("delta_clt with identity g and dg recovers clt result", {
  # For identity g, delta_clt should give the same result as clt
  e <- exponential(rate = 3)
  g <- function(x) x
  dg <- function(x) 1

  result_delta <- delta_clt(e, g, dg)
  result_clt <- clt(e)

  expect_equal(mean(result_delta), mean(result_clt))
  expect_equal(vcov(result_delta), vcov(result_clt), tolerance = 1e-10)
})

test_that("delta_clt errors on non-dist input", {
  expect_error(delta_clt(42, identity, function(x) 1),
               "'base_dist' must be a 'dist' object")
  expect_error(delta_clt("hello", identity, function(x) 1),
               "'base_dist' must be a 'dist' object")
})

test_that("delta_clt handles zero derivative (degenerate case)", {
  # Given: g(x) = constant, so dg(x) = 0
  #   Limiting variance = 0^2 * var = 0
  n <- normal(mu = 3, var = 5)
  g <- function(x) 42
  dg <- function(x) 0

  result <- delta_clt(n, g, dg)

  expect_true(is_normal(result))
  expect_equal(mean(result), 0)
  expect_equal(vcov(result), 0)
})

# =============================================================================
# normal_approx() — Moment-matching normal approximation
# =============================================================================

test_that("normal_approx of gamma distribution has correct moments", {
  # Given: Gamma(shape=3, rate=2), mean = 1.5, var = 0.75
  g <- gamma_dist(shape = 3, rate = 2)

  # When: Computing the normal approximation
  result <- normal_approx(g)

  # Then: normal(1.5, 0.75)
  expect_true(is_normal(result))
  expect_equal(mean(result), 1.5)
  expect_equal(vcov(result), 0.75)
})

test_that("normal_approx of mvn returns equivalent mvn (idempotent)", {
  # Given: A 2D MVN with non-trivial covariance
  sig <- matrix(c(1, 0.5, 0.5, 2), 2, 2)
  m <- mvn(mu = c(1, 2), sigma = sig)

  # When: Computing the normal approximation
  result <- normal_approx(m)

  # Then: The result is an mvn with the same mean and covariance
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(1, 2))
  expect_equal(vcov(result), sig)
})

test_that("normal_approx of normal returns equivalent normal (idempotent)", {
  # Given: A normal distribution
  n <- normal(mu = 7, var = 4)

  # When: Computing the normal approximation
  result <- normal_approx(n)

  # Then: The result has the same parameters
  expect_true(is_normal(result))
  expect_equal(mean(result), 7)
  expect_equal(vcov(result), 4)
})

test_that("normal_approx of exponential has correct moments", {
  # Given: Exponential(rate=5), mean = 0.2, var = 0.04
  e <- exponential(rate = 5)

  # When: Computing the normal approximation
  result <- normal_approx(e)

  # Then: normal(0.2, 0.04)
  expect_true(is_normal(result))
  expect_equal(mean(result), 0.2)
  expect_equal(vcov(result), 0.04)
})

test_that("normal_approx errors on non-dist input", {
  expect_error(normal_approx(42), "'x' must be a 'dist' object")
  expect_error(normal_approx("hello"), "'x' must be a 'dist' object")
  expect_error(normal_approx(list(mu = 0, var = 1)), "'x' must be a 'dist' object")
})

test_that("normal_approx of beta distribution has correct moments", {
  # Given: Beta(shape1=2, shape2=5), mean = 2/7, var = (2*5) / (7^2 * 8)
  b <- beta_dist(shape1 = 2, shape2 = 5)
  expected_mean <- 2 / 7
  expected_var <- (2 * 5) / (7^2 * 8)

  # When: Computing the normal approximation
  result <- normal_approx(b)

  # Then: normal with matching moments
  expect_true(is_normal(result))
  expect_equal(mean(result), expected_mean, tolerance = 1e-10)
  expect_equal(vcov(result), expected_var, tolerance = 1e-10)
})

# =============================================================================
# Integration tests — relationships between limit functions
# =============================================================================

test_that("clt and lln give same mean for degenerate distributions", {
  # LLN gives degenerate at mu; CLT centers at 0.
  # They agree in the sense that both give normal family outputs.
  n <- normal(mu = 10, var = 4)

  lln_result <- lln(n)
  clt_result <- clt(n)

  # LLN has the population mean; CLT is centered at 0

  expect_equal(mean(lln_result), 10)
  expect_equal(mean(clt_result), 0)

  # LLN has zero variance; CLT has population variance
  expect_equal(vcov(lln_result), 0)
  expect_equal(vcov(clt_result), 4)
})

test_that("normal_approx agrees with clt up to location shift", {
  # normal_approx gives N(mu, var); clt gives N(0, var)
  # So normal_approx = location-shifted clt
  e <- exponential(rate = 3)

  approx_result <- normal_approx(e)
  clt_result <- clt(e)

  # Same variance
  expect_equal(vcov(approx_result), vcov(clt_result), tolerance = 1e-10)
  # Means differ by population mean
  expect_equal(mean(approx_result), mean(e))
  expect_equal(mean(clt_result), 0)
})
