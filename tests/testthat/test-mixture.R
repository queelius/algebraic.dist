# Tests for the mixture distribution class

# --- Construction tests ---

test_that("mixture constructor creates valid object with valid args", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 2)
  m <- mixture(list(n1, n2), c(0.3, 0.7))

  expect_s3_class(m, "mixture")
  expect_s3_class(m, "dist")
  expect_length(m$components, 2)
  expect_equal(m$weights, c(0.3, 0.7))
})

test_that("mixture with empty components errors", {
  expect_error(mixture(list(), numeric(0)),
               "'components' must be a non-empty list")
})

test_that("mixture with mismatched weights length errors", {
  n1 <- normal()
  n2 <- normal(mu = 1)
  expect_error(mixture(list(n1, n2), c(0.5, 0.3, 0.2)),
               "'weights' must be a numeric vector of length 2")
})

test_that("mixture with negative weights errors", {
  n1 <- normal()
  n2 <- normal(mu = 1)
  expect_error(mixture(list(n1, n2), c(-0.3, 1.3)),
               "'weights' must be non-negative")
})

test_that("mixture with weights not summing to 1 errors", {
  n1 <- normal()
  n2 <- normal(mu = 1)
  expect_error(mixture(list(n1, n2), c(0.3, 0.3)),
               "'weights' must sum to 1")
})

test_that("mixture with non-dist components errors", {
  expect_error(mixture(list(1, 2), c(0.5, 0.5)),
               "all components must be 'dist' objects")
})

test_that("class hierarchy for continuous univariate mixture", {
  n1 <- normal(mu = 0, var = 1)
  e1 <- exponential(rate = 2)
  m <- mixture(list(n1, e1), c(0.5, 0.5))

  expect_s3_class(m, "mixture")
  expect_s3_class(m, "univariate_dist")
  expect_s3_class(m, "continuous_dist")
  expect_s3_class(m, "dist")
})

test_that("class hierarchy for non-homogeneous components omits sub-types", {
  # A univariate continuous component mixed with a discrete one should

  # not inherit either continuous_dist or discrete_dist
  n1 <- normal(mu = 0, var = 1)
  ed <- empirical_dist(c(1, 2, 3))
  m <- mixture(list(n1, ed), c(0.5, 0.5))

  expect_s3_class(m, "mixture")
  expect_s3_class(m, "univariate_dist")
  expect_s3_class(m, "dist")
  expect_false(inherits(m, "continuous_dist"))
  expect_false(inherits(m, "discrete_dist"))
})

# --- is_mixture ---

test_that("is_mixture works", {
  n1 <- normal()
  n2 <- normal(mu = 3)
  m <- mixture(list(n1, n2), c(0.5, 0.5))

  expect_true(is_mixture(m))
  expect_false(is_mixture(n1))
  expect_false(is_mixture("mixture"))
})

# --- mean ---

test_that("mean() is weighted sum of means", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 10, var = 1)
  m <- mixture(list(n1, n2), c(0.3, 0.7))

  expected_mean <- 0.3 * 0 + 0.7 * 10
  expect_equal(mean(m), expected_mean)
})

test_that("mean() with equal weights is average of component means", {
  e1 <- exponential(rate = 1)
  e2 <- exponential(rate = 0.5)
  m <- mixture(list(e1, e2), c(0.5, 0.5))

  expected_mean <- 0.5 * (1/1) + 0.5 * (1/0.5)
  expect_equal(mean(m), expected_mean)
})

# --- vcov ---

test_that("vcov() uses law of total variance", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 10, var = 4)
  w <- c(0.3, 0.7)
  m <- mixture(list(n1, n2), w)

  comp_means <- c(0, 10)
  comp_vars <- c(1, 4)
  overall_mean <- sum(w * comp_means)
  within_var <- sum(w * comp_vars)
  between_var <- sum(w * (comp_means - overall_mean)^2)
  expected_var <- within_var + between_var

  expect_equal(vcov(m), expected_var)
})

test_that("vcov() of equal normals equals component variance", {
  # Same normal mixed with itself: variance = component variance
  n1 <- normal(mu = 5, var = 3)
  m <- mixture(list(n1, n1), c(0.4, 0.6))

  # between_var = 0 since all means are the same
  # within_var = 3
  expect_equal(vcov(m), 3)
})

# --- density ---

test_that("density() is weighted sum of component densities at several points", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 1)
  w <- c(0.3, 0.7)
  m <- mixture(list(n1, n2), w)

  f <- density(m)
  f1 <- density(n1)
  f2 <- density(n2)

  test_points <- c(-3, -1, 0, 2.5, 5, 8)
  for (t in test_points) {
    expected <- w[1] * f1(t) + w[2] * f2(t)
    expect_equal(f(t), expected, tolerance = 1e-12,
                 label = paste("density at t =", t))
  }
})

test_that("density() handles log argument", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 1)
  m <- mixture(list(n1, n2), c(0.5, 0.5))
  f <- density(m)

  val <- f(2.5)
  log_val <- f(2.5, log = TRUE)
  expect_equal(log_val, log(val), tolerance = 1e-12)
})

test_that("density() is vectorized over t", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 3, var = 1)
  m <- mixture(list(n1, n2), c(0.5, 0.5))
  f <- density(m)

  ts <- c(-1, 0, 1, 2, 3, 4)
  result <- f(ts)
  expect_length(result, length(ts))

  # Check each element matches scalar call
  for (i in seq_along(ts)) {
    expect_equal(result[i], f(ts[i]), tolerance = 1e-12)
  }
})

# --- cdf ---

test_that("cdf() is weighted sum of component CDFs", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 1)
  w <- c(0.3, 0.7)
  m <- mixture(list(n1, n2), w)

  F_mix <- cdf(m)
  F1 <- cdf(n1)
  F2 <- cdf(n2)

  test_points <- c(-3, 0, 2.5, 5, 8)
  for (q in test_points) {
    expected <- w[1] * F1(q) + w[2] * F2(q)
    expect_equal(F_mix(q), expected, tolerance = 1e-12,
                 label = paste("CDF at q =", q))
  }
})

test_that("cdf() is vectorized over q", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 1)
  m <- mixture(list(n1, n2), c(0.5, 0.5))
  F_mix <- cdf(m)

  qs <- c(-2, 0, 2, 4, 6)
  result <- F_mix(qs)
  expect_length(result, length(qs))
})

# --- sampler ---

test_that("sampler() gives correct number of samples", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 1)
  m <- mixture(list(n1, n2), c(0.5, 0.5))

  samp_fn <- sampler(m)
  expect_type(samp_fn, "closure")

  samples <- samp_fn(100)
  expect_length(samples, 100)
  expect_type(samples, "double")
})

test_that("sampler() mean near theoretical mixture mean", {
  set.seed(42)
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 10, var = 1)
  w <- c(0.3, 0.7)
  m <- mixture(list(n1, n2), w)

  samples <- sampler(m)(50000)
  sample_mean <- sum(samples) / length(samples)
  expected_mean <- 0.3 * 0 + 0.7 * 10
  expect_equal(sample_mean, expected_mean, tolerance = 0.1)
})

# --- params ---

test_that("params() concatenates component params and weights", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 2)
  w <- c(0.3, 0.7)
  m <- mixture(list(n1, n2), w)

  p <- params(m)
  expect_true(is.numeric(p))

  # Should contain mu and var from both normals plus the 2 weights
  expect_true("mu" %in% names(p) || "mu1" %in% names(p) || any(grepl("mu", names(p))))
  expect_true(any(grepl("weight", names(p))))
  expect_equal(length(p), 2 + 2 + 2)  # 2 params each + 2 weights
})

# --- nparams ---

test_that("nparams() counts correctly", {
  n1 <- normal(mu = 0, var = 1)
  e1 <- exponential(rate = 2)
  m <- mixture(list(n1, e1), c(0.5, 0.5))

  # normal has 2 params, exponential has 1, plus 2 weights = 5
  expect_equal(nparams(m), 5)
})

# --- dim ---

test_that("dim() returns correct value", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 1)
  m <- mixture(list(n1, n2), c(0.5, 0.5))

  expect_equal(dim(m), 1)
})

# --- sup ---

test_that("sup() spans all component supports", {
  e1 <- exponential(rate = 1)
  n1 <- normal(mu = 0, var = 1)
  m <- mixture(list(e1, n1), c(0.5, 0.5))

  s <- sup(m)
  expect_s3_class(s, "interval")

  # Normal has support (-Inf, Inf), exponential has (0, Inf)
  # So mixture support should be (-Inf, Inf)
  expect_equal(s$infimum(), -Inf)
  expect_equal(s$supremum(), Inf)
})

test_that("sup() with bounded distributions", {
  # Two exponentials: both (0, Inf)
  e1 <- exponential(rate = 1)
  e2 <- exponential(rate = 2)
  m <- mixture(list(e1, e2), c(0.5, 0.5))

  s <- sup(m)
  expect_equal(s$infimum(), 0)
  expect_equal(s$supremum(), Inf)
})

# --- format and print ---

test_that("format() works", {
  n1 <- normal()
  n2 <- normal(mu = 5)
  m <- mixture(list(n1, n2), c(0.5, 0.5))

  fmt <- format(m)
  expect_type(fmt, "character")
  expect_true(grepl("Mixture", fmt))
  expect_true(grepl("2 components", fmt))
})

test_that("format() shows correct number of components", {
  n1 <- normal()
  n2 <- normal(mu = 3)
  n3 <- exponential(rate = 1)
  m <- mixture(list(n1, n2, n3), c(0.2, 0.3, 0.5))

  expect_true(grepl("3 components", format(m)))
})

test_that("print() works", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 5, var = 2)
  m <- mixture(list(n1, n2), c(0.3, 0.7))

  expect_output(print(m), "Mixture distribution")
  expect_output(print(m), "w=0.300")
  expect_output(print(m), "w=0.700")
})

# --- Bimodal density test ---

test_that("mixture of 2 normals has bimodal density", {
  n1 <- normal(mu = -5, var = 1)
  n2 <- normal(mu = 5, var = 1)
  m <- mixture(list(n1, n2), c(0.5, 0.5))

  f <- density(m)

  # Density at the two modes should be higher than at the midpoint
  dens_mode1 <- f(-5)
  dens_mode2 <- f(5)
  dens_midpoint <- f(0)

  expect_gt(dens_mode1, dens_midpoint)
  expect_gt(dens_mode2, dens_midpoint)
})

# --- Mixture of exponentials: CDF between component CDFs ---

test_that("mixture of exponentials: CDF between component CDFs", {
  e1 <- exponential(rate = 1)
  e2 <- exponential(rate = 3)
  m <- mixture(list(e1, e2), c(0.5, 0.5))

  F_mix <- cdf(m)
  F1 <- cdf(e1)
  F2 <- cdf(e2)

  test_points <- c(0.1, 0.5, 1, 2, 5)
  for (t in test_points) {
    val <- F_mix(t)
    lo <- min(F1(t), F2(t))
    hi <- max(F1(t), F2(t))
    expect_true(val >= lo - 1e-12 && val <= hi + 1e-12,
                label = paste("CDF at t =", t, "is between component CDFs"))
  }
})

# --- Edge cases ---

test_that("mixture with a single component behaves like that component", {
  n1 <- normal(mu = 3, var = 2)
  m <- mixture(list(n1), 1)

  expect_equal(mean(m), 3)
  expect_equal(vcov(m), 2)

  f_mix <- density(m)
  f_n1 <- density(n1)
  expect_equal(f_mix(1), f_n1(1), tolerance = 1e-12)
  expect_equal(f_mix(5), f_n1(5), tolerance = 1e-12)
})

test_that("mixture with zero weight on a component ignores it", {
  n1 <- normal(mu = 0, var = 1)
  n2 <- normal(mu = 100, var = 1)
  m <- mixture(list(n1, n2), c(1, 0))

  expect_equal(mean(m), 0)
  expect_equal(vcov(m), 1)
})

test_that("mixture weights tolerance allows tiny numerical error", {
  n1 <- normal()
  n2 <- normal(mu = 1)
  # Weights that sum to 1 within tolerance
  w <- c(1/3, 2/3)
  expect_no_error(mixture(list(n1, n2), w))
})

# --- Multivariate mixture tests ---

test_that("multivariate mixture sampler returns correct shape", {
  set.seed(42)
  m1 <- mvn(mu = c(0, 0), sigma = diag(2))
  m2 <- mvn(mu = c(5, 5), sigma = diag(2))
  mix <- mixture(list(m1, m2), c(0.5, 0.5))

  samp_fn <- sampler(mix)
  samples <- samp_fn(100)
  expect_true(is.matrix(samples))
  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 2)
})

test_that("multivariate mixture vcov returns correct covariance matrix", {
  m1 <- mvn(mu = c(0, 0), sigma = diag(2))
  m2 <- mvn(mu = c(4, 4), sigma = diag(2))
  mix <- mixture(list(m1, m2), c(0.5, 0.5))

  V <- vcov(mix)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 2)
  expect_equal(ncol(V), 2)
  # Within-component variance: 0.5 * I + 0.5 * I = I
  # Between-component variance: 0.5 * outer(c(-2,-2), c(-2,-2)) +
  #                              0.5 * outer(c(2,2), c(2,2)) = 4 * ones
  # Total diagonal should be 1 + 4 = 5
  expect_equal(V[1, 1], 5, tolerance = 1e-10)
  expect_equal(V[2, 2], 5, tolerance = 1e-10)
  expect_equal(V[1, 2], 4, tolerance = 1e-10)
})

test_that("multivariate mixture inherits multivariate_dist class", {
  m1 <- mvn(mu = c(0, 0))
  m2 <- mvn(mu = c(1, 1))
  mix <- mixture(list(m1, m2), c(0.5, 0.5))

  expect_s3_class(mix, "multivariate_dist")
})
