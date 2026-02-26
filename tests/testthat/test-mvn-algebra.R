# Tests for MVN algebra: conditional, affine_transform

# =============================================================================
# conditional.mvn — Schur complement
# =============================================================================

test_that("conditional.mvn with given_indices returns correct type", {
  m <- mvn(mu = c(1, 2, 3), sigma = diag(3))
  # Condition on variable 3 = 5
  result <- conditional(m, given_indices = 3, given_values = 5)
  expect_true(is_mvn(result))
  expect_equal(dim(result), 2)
})

test_that("conditional.mvn independent components: conditioning doesn't affect mean", {
  # Independent components: Sigma = I
  m <- mvn(mu = c(10, 20, 30), sigma = diag(3))
  # Condition on X3 = 50
  result <- conditional(m, given_indices = 3, given_values = 50)
  # Since independent, X1|X3 should still have mean 10, X2|X3 mean 20
  expect_equal(mean(result), c(10, 20))
  expect_equal(vcov(result), diag(2))
})

test_that("conditional.mvn correlated: Schur complement correct", {
  # X = (X1, X2) with cov(X1, X2) = 0.5
  sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  m <- mvn(mu = c(0, 0), sigma = sigma)

  # Condition on X2 = 2
  result <- conditional(m, given_indices = 2, given_values = 2)

  # Expected: mu_cond = 0 + 0.5 * 1 * (2 - 0) = 1
  # sig_cond = 1 - 0.5 * 1 * 0.5 = 0.75
  expect_true(is_normal(result))
  expect_equal(mean(result), 1)
  expect_equal(vcov(result), 0.75)
})

test_that("conditional.mvn 4D: condition on 2 variables", {
  mu <- c(1, 2, 3, 4)
  sigma <- diag(4)
  sigma[1, 3] <- sigma[3, 1] <- 0.8
  sigma[2, 4] <- sigma[4, 2] <- 0.6
  m <- mvn(mu = mu, sigma = sigma)

  # Condition on X3 = 5, X4 = 6
  result <- conditional(m, given_indices = c(3, 4), given_values = c(5, 6))
  expect_true(is_mvn(result))
  expect_equal(dim(result), 2)

  # mu_cond for X1: 1 + 0.8 * 1 * (5 - 3) = 1 + 1.6 = 2.6
  # mu_cond for X2: 2 + 0.6 * 1 * (6 - 4) = 2 + 1.2 = 3.2
  expect_equal(mean(result)[1], 2.6, tolerance = 1e-10)
  expect_equal(mean(result)[2], 3.2, tolerance = 1e-10)
})

test_that("conditional.mvn condition on all-but-one returns normal", {
  m <- mvn(mu = c(0, 0), sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2))
  result <- conditional(m, given_indices = 1, given_values = 1)
  expect_true(is_normal(result))
  expect_equal(mean(result), 0.5)  # 0 + 0.5 * 1 * (1 - 0)
  expect_equal(vcov(result), 0.75) # 1 - 0.5^2
})

test_that("conditional.mvn errors: condition on all variables", {
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  expect_error(
    conditional(m, given_indices = c(1, 2), given_values = c(3, 4)),
    "cannot condition on all variables"
  )
})

test_that("conditional.mvn errors: mismatched indices/values", {
  m <- mvn(mu = c(1, 2, 3), sigma = diag(3))
  expect_error(
    conditional(m, given_indices = c(1, 2), given_values = 3),
    "same length"
  )
})

test_that("conditional.mvn errors: out-of-range indices", {
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  expect_error(
    conditional(m, given_indices = 3, given_values = 1),
    "given_indices"
  )
})

test_that("conditional.mvn P fallback works", {
  m <- mvn(mu = c(0, 0), sigma = diag(2))
  set.seed(42)
  result <- conditional(m, P = function(x) x[1] > 0)
  expect_s3_class(result, "empirical_dist")
})

test_that("conditional.mvn errors: no P or given_indices", {
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  expect_error(conditional(m), "must provide either")
})

# ---- MC validation of Schur complement ----

test_that("conditional.mvn: MC validation of closed form", {
  sigma <- matrix(c(4, 2, 2, 3), 2, 2)
  m <- mvn(mu = c(5, 10), sigma = sigma)

  # Closed form: X1 | X2 = 12
  result <- conditional(m, given_indices = 2, given_values = 12)

  # MC validation: sample from joint, filter X2 ≈ 12, check X1 distribution
  set.seed(42)
  samp <- sampler(m)(100000)
  # Select samples where X2 is close to 12
  close <- abs(samp[, 2] - 12) < 0.5
  x1_given <- samp[close, 1]

  # Closed-form mean should be close to empirical mean
  expect_equal(mean(result), mean(x1_given), tolerance = 0.3)
})


# =============================================================================
# affine_transform
# =============================================================================

test_that("affine_transform: identity matrix on mvn returns same distribution", {
  m <- mvn(mu = c(1, 2), sigma = matrix(c(4, 1, 1, 3), 2, 2))
  result <- affine_transform(m, A = diag(2))
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(1, 2))
  expect_equal(vcov(result), matrix(c(4, 1, 1, 3), 2, 2))
})

test_that("affine_transform: scale matrix", {
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  A <- matrix(c(2, 0, 0, 3), 2, 2)
  result <- affine_transform(m, A = A)
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(2, 6))
  expect_equal(vcov(result), matrix(c(4, 0, 0, 9), 2, 2))
})

test_that("affine_transform: with offset b", {
  m <- mvn(mu = c(0, 0), sigma = diag(2))
  result <- affine_transform(m, A = diag(2), b = c(5, 10))
  expect_true(is_mvn(result))
  expect_equal(mean(result), c(5, 10))
  expect_equal(vcov(result), diag(2))
})

test_that("affine_transform: projection (2D -> 1D)", {
  m <- mvn(mu = c(3, 7), sigma = matrix(c(4, 1, 1, 2), 2, 2))
  # Project onto first component
  A <- matrix(c(1, 0), nrow = 1)
  result <- affine_transform(m, A = A)
  expect_true(is_normal(result))
  expect_equal(mean(result), 3)
  expect_equal(vcov(result), 4)
})

test_that("affine_transform: sum of components (1xd matrix)", {
  m <- mvn(mu = c(1, 2), sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2))
  # A = [1, 1] -> X1 + X2
  A <- matrix(c(1, 1), nrow = 1)
  result <- affine_transform(m, A = A)
  expect_true(is_normal(result))
  expect_equal(mean(result), 3)
  # var(X1 + X2) = 1 + 2*0.5 + 1 = 3
  expect_equal(vcov(result), 3)
})

test_that("affine_transform: works on univariate normal", {
  n <- normal(mu = 5, var = 4)
  # 2*X + 3
  result <- affine_transform(n, A = 2, b = 3)
  expect_true(is_normal(result))
  expect_equal(mean(result), 13)   # 2*5 + 3
  expect_equal(vcov(result), 16)   # 2^2 * 4
})

test_that("affine_transform errors: wrong dimensions", {
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  expect_error(affine_transform(m, A = diag(3)), "ncol.*must equal dim")
})

test_that("affine_transform errors: wrong b length", {
  m <- mvn(mu = c(1, 2), sigma = diag(2))
  expect_error(affine_transform(m, A = diag(2), b = c(1, 2, 3)), "length.*must equal nrow")
})

test_that("affine_transform errors: non-dist input", {
  expect_error(affine_transform(42, A = 1), "'x' must be a 'dist' object")
})

test_that("affine_transform errors: unsupported dist type", {
  expect_error(affine_transform(exponential(1), A = 1), "'x' must be a 'normal' or 'mvn'")
})
