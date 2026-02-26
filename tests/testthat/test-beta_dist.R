# Tests for the beta distribution class

# --- Construction ---

test_that("beta_dist constructor creates valid object with correct params", {
  b <- beta_dist(shape1 = 2, shape2 = 5)

  expect_s3_class(b, "beta_dist")
  expect_s3_class(b, "univariate_dist")
  expect_s3_class(b, "continuous_dist")
  expect_s3_class(b, "dist")
  expect_equal(b$shape1, 2)
  expect_equal(b$shape2, 5)
})

test_that("beta_dist constructor rejects invalid shape1", {
  expect_error(beta_dist(shape1 = -1, shape2 = 1), "'shape1' must be a positive scalar")
  expect_error(beta_dist(shape1 = 0, shape2 = 1), "'shape1' must be a positive scalar")
  expect_error(beta_dist(shape1 = "a", shape2 = 1), "'shape1' must be a positive scalar")
  expect_error(beta_dist(shape1 = c(1, 2), shape2 = 1), "'shape1' must be a positive scalar")
  expect_error(beta_dist(shape1 = NA_real_, shape2 = 1), "'shape1' must be a positive scalar")
})

test_that("beta_dist constructor rejects invalid shape2", {
  expect_error(beta_dist(shape1 = 1, shape2 = -1), "'shape2' must be a positive scalar")
  expect_error(beta_dist(shape1 = 1, shape2 = 0), "'shape2' must be a positive scalar")
  expect_error(beta_dist(shape1 = 1, shape2 = "b"), "'shape2' must be a positive scalar")
  expect_error(beta_dist(shape1 = 1, shape2 = c(1, 2)), "'shape2' must be a positive scalar")
  expect_error(beta_dist(shape1 = 1, shape2 = NA_real_), "'shape2' must be a positive scalar")
})

# --- Class hierarchy ---

test_that("beta_dist has correct class hierarchy", {
  b <- beta_dist(shape1 = 2, shape2 = 3)
  classes <- class(b)
  expect_equal(classes, c("beta_dist", "univariate_dist",
                          "continuous_dist", "dist"))
})

# --- is_beta_dist ---

test_that("is_beta_dist identifies beta_dist objects correctly", {
  b <- beta_dist(shape1 = 2, shape2 = 3)

  expect_true(is_beta_dist(b))
  expect_false(is_beta_dist(list(shape1 = 2, shape2 = 3)))
  expect_false(is_beta_dist(normal()))
  expect_false(is_beta_dist(uniform_dist()))
})

# --- params ---

test_that("params.beta_dist returns named vector of parameters", {
  b <- beta_dist(shape1 = 2.5, shape2 = 0.7)
  p <- params(b)

  expect_named(p, c("shape1", "shape2"))
  expect_equal(p["shape1"], c(shape1 = 2.5))
  expect_equal(p["shape2"], c(shape2 = 0.7))
})

# --- mean ---

test_that("mean.beta_dist returns shape1/(shape1+shape2)", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  expect_equal(mean(b), 2 / 7)

  b2 <- beta_dist(shape1 = 1, shape2 = 1)
  expect_equal(mean(b2), 0.5)

  b3 <- beta_dist(shape1 = 3, shape2 = 3)
  expect_equal(mean(b3), 0.5)

  b4 <- beta_dist(shape1 = 0.5, shape2 = 0.5)
  expect_equal(mean(b4), 0.5)
})

# --- vcov ---

test_that("vcov.beta_dist returns shape1*shape2/((a+b)^2*(a+b+1))", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  expected <- 2 * 5 / ((2 + 5)^2 * (2 + 5 + 1))
  expect_equal(vcov(b), expected)

  b2 <- beta_dist(shape1 = 1, shape2 = 1)
  expected2 <- 1 * 1 / ((1 + 1)^2 * (1 + 1 + 1))
  expect_equal(vcov(b2), expected2)
  expect_equal(vcov(b2), 1 / 12)  # Uniform(0,1) variance
})

# --- dim ---

test_that("dim.beta_dist returns 1 for univariate distribution", {
  b <- beta_dist(shape1 = 2, shape2 = 3)
  expect_equal(dim(b), 1)
})

# --- format / print ---

test_that("format.beta_dist returns correct string", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  expect_equal(format(b), "Beta distribution (shape1 = 2, shape2 = 5)")
})

test_that("print.beta_dist outputs to console", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  expect_output(print(b), "Beta distribution")
  expect_output(print(b), "shape1 = 2")
  expect_output(print(b), "shape2 = 5")
})

test_that("print.beta_dist returns object invisibly", {
  b <- beta_dist(shape1 = 1, shape2 = 1)
  out <- capture.output(ret <- print(b))
  expect_identical(ret, b)
})

# --- sampler ---

test_that("sampler.beta_dist returns a function that generates samples", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  samp_fn <- sampler(b)

  expect_type(samp_fn, "closure")

  samples <- samp_fn(100)
  expect_length(samples, 100)
})

test_that("sampler.beta_dist produces samples within (0, 1)", {
  set.seed(42)
  b <- beta_dist(shape1 = 2, shape2 = 5)
  samples <- sampler(b)(10000)

  expect_true(all(samples > 0))
  expect_true(all(samples < 1))
})

test_that("sampler.beta_dist produces samples with approximately correct mean", {
  set.seed(42)
  b <- beta_dist(shape1 = 2, shape2 = 5)
  samples <- sampler(b)(10000)

  expected_mean <- 2 / (2 + 5)
  sample_mean <- sum(samples) / length(samples)
  expect_equal(sample_mean, expected_mean, tolerance = 0.05)
})

# --- density ---

test_that("density.beta_dist returns correct probability density", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  pdf <- density(b)

  # Compare with dbeta at known points
  expect_equal(pdf(0.3), dbeta(0.3, shape1 = 2, shape2 = 5), tolerance = 1e-10)
  expect_equal(pdf(0.5), dbeta(0.5, shape1 = 2, shape2 = 5), tolerance = 1e-10)
  expect_equal(pdf(0.9), dbeta(0.9, shape1 = 2, shape2 = 5), tolerance = 1e-10)
})

test_that("density.beta_dist handles log argument correctly", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  pdf <- density(b)

  expect_equal(pdf(0.3, log = TRUE),
               dbeta(0.3, shape1 = 2, shape2 = 5, log = TRUE),
               tolerance = 1e-10)
})

test_that("density.beta_dist returns zero outside (0, 1)", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  pdf <- density(b)

  expect_equal(pdf(-0.1), 0)
  expect_equal(pdf(1.1), 0)
})

# --- cdf ---

test_that("cdf.beta_dist returns correct cumulative distribution", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  cdf_fn <- cdf(b)

  expect_equal(cdf_fn(0), pbeta(0, shape1 = 2, shape2 = 5), tolerance = 1e-10)
  expect_equal(cdf_fn(0.5), pbeta(0.5, shape1 = 2, shape2 = 5), tolerance = 1e-10)
  expect_equal(cdf_fn(1), pbeta(1, shape1 = 2, shape2 = 5), tolerance = 1e-10)
  expect_equal(cdf_fn(-0.1), 0, tolerance = 1e-10)
  expect_equal(cdf_fn(1.1), 1, tolerance = 1e-10)
})

test_that("cdf.beta_dist handles log.p argument correctly", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  cdf_fn <- cdf(b)

  expect_equal(cdf_fn(0.5, log.p = TRUE),
               pbeta(0.5, shape1 = 2, shape2 = 5, log.p = TRUE),
               tolerance = 1e-10)
})

# --- inv_cdf ---

test_that("inv_cdf.beta_dist returns correct quantiles", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  qf <- inv_cdf(b)

  expect_equal(qf(0.5), qbeta(0.5, shape1 = 2, shape2 = 5), tolerance = 1e-10)
  expect_equal(qf(0.95), qbeta(0.95, shape1 = 2, shape2 = 5), tolerance = 1e-10)
  expect_equal(qf(0.1), qbeta(0.1, shape1 = 2, shape2 = 5), tolerance = 1e-10)
})

test_that("inv_cdf.beta_dist round-trips with cdf", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  cdf_fn <- cdf(b)
  qf <- inv_cdf(b)

  # cdf(inv_cdf(p)) == p
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (p in probs) {
    expect_equal(cdf_fn(qf(p)), p, tolerance = 1e-10)
  }

  # inv_cdf(cdf(x)) == x for x in (0, 1)
  vals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  for (v in vals) {
    expect_equal(qf(cdf_fn(v)), v, tolerance = 1e-10)
  }
})

# --- sup ---

test_that("sup.beta_dist returns open interval (0, 1)", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  s <- sup(b)

  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), 0)
  expect_equal(s$supremum(), 1)
  expect_false(s$lower_closed)
  expect_false(s$upper_closed)
})

test_that("sup.beta_dist support contains interior but not boundary points", {
  b <- beta_dist(shape1 = 2, shape2 = 5)
  s <- sup(b)

  expect_true(s$has(0.5))
  expect_true(s$has(0.01))
  expect_true(s$has(0.99))
  expect_false(s$has(0))
  expect_false(s$has(1))
  expect_false(s$has(-0.1))
  expect_false(s$has(1.1))
})

# --- Cross-validation: Beta(1,1) == Uniform(0,1) ---

test_that("beta_dist(1, 1) density matches uniform_dist(0, 1) density", {
  b <- beta_dist(shape1 = 1, shape2 = 1)
  u <- uniform_dist(min = 0, max = 1)
  b_pdf <- density(b)
  u_pdf <- density(u)

  # Beta(1,1) is Uniform(0,1): density should be 1 on (0,1)
  test_points <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  for (t in test_points) {
    expect_equal(b_pdf(t), u_pdf(t), tolerance = 1e-10,
                 label = paste("density at t =", t))
    expect_equal(b_pdf(t), 1, tolerance = 1e-10,
                 label = paste("density = 1 at t =", t))
  }
})

test_that("beta_dist(1, 1) CDF matches uniform_dist(0, 1) CDF", {
  b <- beta_dist(shape1 = 1, shape2 = 1)
  u <- uniform_dist(min = 0, max = 1)
  b_cdf <- cdf(b)
  u_cdf <- cdf(u)

  test_points <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  for (t in test_points) {
    expect_equal(b_cdf(t), u_cdf(t), tolerance = 1e-10,
                 label = paste("CDF at t =", t))
  }
})

test_that("beta_dist(1, 1) mean and variance match Uniform(0,1)", {
  b <- beta_dist(shape1 = 1, shape2 = 1)
  u <- uniform_dist(min = 0, max = 1)

  expect_equal(mean(b), mean(u), tolerance = 1e-10)
  expect_equal(vcov(b), vcov(u), tolerance = 1e-10)
})
