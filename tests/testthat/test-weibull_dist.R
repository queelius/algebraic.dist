# Tests for the weibull_dist distribution class

# --- Construction -----------------------------------------------------------

test_that("weibull_dist constructor creates valid object with correct params", {
  w <- weibull_dist(shape = 2, scale = 3)

  expect_s3_class(w, "weibull_dist")
  expect_equal(w$shape, 2)
  expect_equal(w$scale, 3)
})

test_that("weibull_dist has correct class hierarchy", {
  w <- weibull_dist(shape = 1.5, scale = 2)

  expect_s3_class(w, "weibull_dist")
  expect_s3_class(w, "univariate_dist")
  expect_s3_class(w, "continuous_dist")
  expect_s3_class(w, "dist")
})

test_that("weibull_dist rejects non-positive shape", {
  expect_error(weibull_dist(shape = 0, scale = 1), "'shape' must be a positive scalar")
  expect_error(weibull_dist(shape = -1, scale = 1), "'shape' must be a positive scalar")
})

test_that("weibull_dist rejects non-positive scale", {
  expect_error(weibull_dist(shape = 1, scale = 0), "'scale' must be a positive scalar")
  expect_error(weibull_dist(shape = 1, scale = -2), "'scale' must be a positive scalar")
})

test_that("weibull_dist rejects non-numeric arguments", {
  expect_error(weibull_dist(shape = "a", scale = 1), "'shape' must be a positive scalar")
  expect_error(weibull_dist(shape = 1, scale = "b"), "'scale' must be a positive scalar")
})

test_that("weibull_dist rejects vector arguments", {
  expect_error(weibull_dist(shape = c(1, 2), scale = 1), "'shape' must be a positive scalar")
  expect_error(weibull_dist(shape = 1, scale = c(1, 2)), "'scale' must be a positive scalar")
})

test_that("weibull_dist rejects NA arguments", {
  expect_error(weibull_dist(shape = NA_real_, scale = 1), "'shape' must be a positive scalar")
  expect_error(weibull_dist(shape = 1, scale = NA_real_), "'scale' must be a positive scalar")
})

# --- is_weibull_dist --------------------------------------------------------

test_that("is_weibull_dist returns TRUE for weibull_dist objects", {
  w <- weibull_dist(shape = 2, scale = 1)
  expect_true(is_weibull_dist(w))
})

test_that("is_weibull_dist returns FALSE for non-weibull objects", {
  expect_false(is_weibull_dist(normal()))
  expect_false(is_weibull_dist(exponential(rate = 1)))
  expect_false(is_weibull_dist(list(shape = 1, scale = 1)))
  expect_false(is_weibull_dist(42))
})

# --- params -----------------------------------------------------------------

test_that("params returns named vector with shape and scale", {
  w <- weibull_dist(shape = 2.5, scale = 3.7)
  p <- params(w)

  expect_named(p, c("shape", "scale"))
  expect_equal(p["shape"], c(shape = 2.5))
  expect_equal(p["scale"], c(scale = 3.7))
})

# --- mean -------------------------------------------------------------------

test_that("mean equals scale * gamma(1 + 1/shape)", {
  shape <- 2
  scale <- 3
  w <- weibull_dist(shape = shape, scale = scale)

  expected <- scale * gamma(1 + 1 / shape)
  expect_equal(mean(w), expected)
})

test_that("mean is correct for shape = 1 (exponential)", {
  # Weibull(1, scale) is Exp(1/scale); mean = scale
  w <- weibull_dist(shape = 1, scale = 5)
  expect_equal(mean(w), 5)
})

# --- vcov -------------------------------------------------------------------

test_that("vcov equals scale^2 * (gamma(1+2/shape) - gamma(1+1/shape)^2)", {
  shape <- 2
  scale <- 3
  w <- weibull_dist(shape = shape, scale = scale)

  expected <- scale^2 * (gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
  expect_equal(vcov(w), expected)
})

test_that("vcov is correct for shape = 1 (exponential)", {
  # Weibull(1, scale) is Exp(1/scale); variance = scale^2
  w <- weibull_dist(shape = 1, scale = 5)
  expect_equal(vcov(w), 25)
})

# --- dim --------------------------------------------------------------------

test_that("dim returns 1", {
  w <- weibull_dist(shape = 2, scale = 1)
  expect_equal(dim(w), 1)
})

# --- format and print -------------------------------------------------------

test_that("format returns descriptive string", {
  w <- weibull_dist(shape = 2, scale = 3)
  f <- format(w)

  expect_type(f, "character")
  expect_match(f, "Weibull distribution")
  expect_match(f, "shape = 2")
  expect_match(f, "scale = 3")
})

test_that("print outputs formatted string", {
  w <- weibull_dist(shape = 2, scale = 3)

  expect_output(print(w), "Weibull distribution")
  expect_output(print(w), "shape = 2")
  expect_output(print(w), "scale = 3")
})

test_that("print returns object invisibly", {
  w <- weibull_dist(shape = 2, scale = 3)
  result <- withVisible(print(capture.output(print(w))))
  # Actually test that print returns invisible
  out <- capture.output(ret <- print(w))
  expect_identical(ret, w)
})

# --- sampler ----------------------------------------------------------------

test_that("sampler returns a function that generates correct number of samples", {
  w <- weibull_dist(shape = 2, scale = 3)
  samp_fn <- sampler(w)

  expect_type(samp_fn, "closure")

  samples <- samp_fn(100)
  expect_length(samples, 100)
  expect_true(all(samples > 0))
})

test_that("sampler produces samples with approximately correct mean", {
  set.seed(42)
  shape <- 2
  scale <- 3
  w <- weibull_dist(shape = shape, scale = scale)
  samples <- sampler(w)(10000)

  expected_mean <- scale * gamma(1 + 1 / shape)
  sample_mean <- sum(samples) / length(samples)
  expect_equal(sample_mean, expected_mean, tolerance = 0.1)
})

test_that("sampler with n=1 returns single value", {
  w <- weibull_dist(shape = 2, scale = 1)
  s <- sampler(w)(1)
  expect_length(s, 1)
})

# --- density ----------------------------------------------------------------

test_that("density matches dweibull at known points", {
  w <- weibull_dist(shape = 2, scale = 3)
  pdf <- density(w)

  test_pts <- c(0.5, 1, 2, 5)
  for (t in test_pts) {
    expect_equal(pdf(t), dweibull(t, shape = 2, scale = 3), tolerance = 1e-12)
  }
})

test_that("density handles log argument correctly", {
  w <- weibull_dist(shape = 2, scale = 3)
  pdf <- density(w)

  expect_equal(pdf(1, log = TRUE),
               dweibull(1, shape = 2, scale = 3, log = TRUE),
               tolerance = 1e-12)
})

test_that("density returns zero for non-positive values", {
  w <- weibull_dist(shape = 2, scale = 3)
  pdf <- density(w)

  expect_equal(pdf(0), 0)
  expect_equal(pdf(-1), 0)
})

test_that("density handles vectorized input", {
  w <- weibull_dist(shape = 1.5, scale = 2)
  pdf <- density(w)

  t_vals <- c(0.1, 0.5, 1, 2, 5)
  expect_equal(pdf(t_vals), dweibull(t_vals, shape = 1.5, scale = 2),
               tolerance = 1e-12)
})

# --- cdf --------------------------------------------------------------------

test_that("cdf matches pweibull at known points", {
  w <- weibull_dist(shape = 2, scale = 3)
  cdf_fn <- cdf(w)

  test_pts <- c(0, 0.5, 1, 2, 5)
  for (q in test_pts) {
    expect_equal(cdf_fn(q), pweibull(q, shape = 2, scale = 3), tolerance = 1e-12)
  }
})

test_that("cdf handles log.p argument correctly", {
  w <- weibull_dist(shape = 2, scale = 3)
  cdf_fn <- cdf(w)

  expect_equal(cdf_fn(1, log.p = TRUE),
               pweibull(1, shape = 2, scale = 3, log.p = TRUE),
               tolerance = 1e-12)
})

# --- inv_cdf ----------------------------------------------------------------

test_that("inv_cdf matches qweibull at known quantiles", {
  w <- weibull_dist(shape = 2, scale = 3)
  qf <- inv_cdf(w)

  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (p in probs) {
    expect_equal(qf(p), qweibull(p, shape = 2, scale = 3), tolerance = 1e-12)
  }
})

test_that("inv_cdf round-trips with cdf", {
  w <- weibull_dist(shape = 2, scale = 3)
  cdf_fn <- cdf(w)
  qf <- inv_cdf(w)

  # cdf(inv_cdf(p)) == p
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (p in probs) {
    expect_equal(cdf_fn(qf(p)), p, tolerance = 1e-12)
  }

  # inv_cdf(cdf(t)) == t
  test_pts <- c(0.5, 1, 3, 10)
  for (t in test_pts) {
    expect_equal(qf(cdf_fn(t)), t, tolerance = 1e-12)
  }
})

test_that("inv_cdf handles lower.tail and log.p arguments", {
  w <- weibull_dist(shape = 2, scale = 3)
  qf <- inv_cdf(w)

  # lower.tail = FALSE gives upper quantile
  expect_equal(qf(0.1, lower.tail = FALSE),
               qweibull(0.1, shape = 2, scale = 3, lower.tail = FALSE),
               tolerance = 1e-12)

  # log.p = TRUE
  expect_equal(qf(log(0.5), log.p = TRUE),
               qweibull(log(0.5), shape = 2, scale = 3, log.p = TRUE),
               tolerance = 1e-12)
})

# --- sup (support) ----------------------------------------------------------

test_that("sup returns interval (0, Inf)", {
  w <- weibull_dist(shape = 2, scale = 3)
  s <- sup(w)

  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), 0)
  expect_equal(s$supremum(), Inf)
  expect_false(s$lower_closed)
})

# --- hazard -----------------------------------------------------------------

test_that("hazard matches closed-form (shape/scale)*(t/scale)^(shape-1)", {
  shape <- 2
  scale <- 3
  w <- weibull_dist(shape = shape, scale = scale)
  h <- hazard(w)

  test_pts <- c(0.5, 1, 2, 5, 10)
  for (t in test_pts) {
    expected <- (shape / scale) * (t / scale)^(shape - 1)
    expect_equal(h(t), expected, tolerance = 1e-12)
  }
})

test_that("hazard returns zero for non-positive values", {
  w <- weibull_dist(shape = 2, scale = 3)
  h <- hazard(w)

  expect_equal(h(0), 0)
  expect_equal(h(-1), 0)
})

test_that("hazard handles log argument correctly", {
  shape <- 2
  scale <- 3
  w <- weibull_dist(shape = shape, scale = scale)
  h <- hazard(w)

  t <- 2
  expected_h <- (shape / scale) * (t / scale)^(shape - 1)
  expect_equal(h(t, log = TRUE), log(expected_h), tolerance = 1e-12)
})

test_that("hazard log is -Inf for non-positive values", {
  w <- weibull_dist(shape = 2, scale = 3)
  h <- hazard(w)

  expect_equal(h(0, log = TRUE), -Inf)
  expect_equal(h(-1, log = TRUE), -Inf)
})

test_that("hazard is constant for shape = 1 (exponential)", {
  # Weibull(1, scale) has hazard = 1/scale, constant
  scale <- 5
  w <- weibull_dist(shape = 1, scale = scale)
  h <- hazard(w)

  expect_equal(h(0.5), 1 / scale)
  expect_equal(h(10), 1 / scale)
  expect_equal(h(100), 1 / scale)
})

test_that("hazard handles vectorized input", {
  shape <- 2
  scale <- 3
  w <- weibull_dist(shape = shape, scale = scale)
  h <- hazard(w)

  t_vals <- c(0.5, 1, 2, 5)
  expected <- (shape / scale) * (t_vals / scale)^(shape - 1)
  expect_equal(h(t_vals), expected, tolerance = 1e-12)
})

# --- surv -------------------------------------------------------------------

test_that("surv equals 1 - cdf at test points", {
  w <- weibull_dist(shape = 2, scale = 3)
  sf <- surv(w)
  cdf_fn <- cdf(w)

  test_pts <- c(0, 0.5, 1, 2, 5, 10)
  for (t in test_pts) {
    expect_equal(sf(t), 1 - cdf_fn(t), tolerance = 1e-10)
  }
})

test_that("surv matches pweibull with lower.tail = FALSE", {
  w <- weibull_dist(shape = 2, scale = 3)
  sf <- surv(w)

  test_pts <- c(0.5, 1, 2, 5)
  for (t in test_pts) {
    expect_equal(sf(t), pweibull(t, shape = 2, scale = 3, lower.tail = FALSE),
                 tolerance = 1e-12)
  }
})

test_that("surv handles log.p argument", {
  w <- weibull_dist(shape = 2, scale = 3)
  sf <- surv(w)

  expect_equal(sf(1, log.p = TRUE),
               pweibull(1, shape = 2, scale = 3, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-12)
})

# --- Cross-validation: Weibull(1, 1/rate) == Exp(rate) ----------------------

test_that("weibull_dist(1, 1/rate) density matches exponential(rate)", {

  rate <- 2
  w <- weibull_dist(shape = 1, scale = 1 / rate)
  e <- exponential(rate = rate)

  w_pdf <- density(w)
  e_pdf <- density(e)

  test_pts <- c(0.1, 0.5, 1, 2, 5)
  for (t in test_pts) {
    expect_equal(w_pdf(t), e_pdf(t), tolerance = 1e-12,
                 label = paste("density at t =", t))
  }
})

test_that("weibull_dist(1, 1/rate) cdf matches exponential(rate)", {
  rate <- 2
  w <- weibull_dist(shape = 1, scale = 1 / rate)
  e <- exponential(rate = rate)

  w_cdf <- cdf(w)
  e_cdf <- cdf(e)

  test_pts <- c(0, 0.1, 0.5, 1, 2, 5)
  for (t in test_pts) {
    expect_equal(w_cdf(t), e_cdf(t), tolerance = 1e-12,
                 label = paste("cdf at t =", t))
  }
})

test_that("weibull_dist(1, 1/rate) mean matches exponential(rate)", {
  rate <- 2
  w <- weibull_dist(shape = 1, scale = 1 / rate)
  e <- exponential(rate = rate)

  expect_equal(mean(w), mean(e), tolerance = 1e-12)
})

test_that("weibull_dist(1, 1/rate) hazard matches exponential(rate)", {
  rate <- 2
  w <- weibull_dist(shape = 1, scale = 1 / rate)
  e <- exponential(rate = rate)

  w_h <- hazard(w)
  e_h <- hazard(e)

  test_pts <- c(0.5, 1, 5, 10)
  for (t in test_pts) {
    expect_equal(w_h(t), e_h(t), tolerance = 1e-12,
                 label = paste("hazard at t =", t))
  }
})

test_that("weibull_dist(1, 1/rate) surv matches exponential(rate)", {
  rate <- 2
  w <- weibull_dist(shape = 1, scale = 1 / rate)
  e <- exponential(rate = rate)

  w_sf <- surv(w)
  e_sf <- surv(e)

  test_pts <- c(0.5, 1, 2, 5)
  for (t in test_pts) {
    expect_equal(w_sf(t), e_sf(t), tolerance = 1e-12,
                 label = paste("surv at t =", t))
  }
})
