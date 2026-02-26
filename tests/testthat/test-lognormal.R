# Tests for the lognormal distribution class

# --- Construction ---

test_that("lognormal constructor creates valid object with default params", {
  d <- lognormal()

  expect_s3_class(d, "lognormal")
  expect_s3_class(d, "univariate_dist")
  expect_s3_class(d, "continuous_dist")
  expect_s3_class(d, "dist")
  expect_equal(d$meanlog, 0)
  expect_equal(d$sdlog, 1)
})

test_that("lognormal constructor creates valid object with custom params", {
  d <- lognormal(meanlog = 2, sdlog = 0.5)

  expect_s3_class(d, "lognormal")
  expect_equal(d$meanlog, 2)
  expect_equal(d$sdlog, 0.5)
})

test_that("lognormal constructor accepts negative meanlog", {
  d <- lognormal(meanlog = -3, sdlog = 1)
  expect_equal(d$meanlog, -3)
})

test_that("lognormal constructor rejects invalid sdlog", {
  expect_error(lognormal(sdlog = 0), "'sdlog' must be a positive scalar")
  expect_error(lognormal(sdlog = -1), "'sdlog' must be a positive scalar")
  expect_error(lognormal(sdlog = "a"), "'sdlog' must be a positive scalar")
  expect_error(lognormal(sdlog = c(1, 2)), "'sdlog' must be a positive scalar")
  expect_error(lognormal(sdlog = NA_real_), "'sdlog' must be a positive scalar")
})

test_that("lognormal constructor rejects non-numeric meanlog", {
  expect_error(lognormal(meanlog = "a"), "'meanlog' must be a numeric scalar")
  expect_error(lognormal(meanlog = c(1, 2)), "'meanlog' must be a numeric scalar")
  expect_error(lognormal(meanlog = NA_real_), "'meanlog' must be a numeric scalar")
})

# --- is_lognormal ---

test_that("is_lognormal identifies lognormal objects correctly", {
  d <- lognormal()

  expect_true(is_lognormal(d))
  expect_false(is_lognormal(list(meanlog = 0, sdlog = 1)))
  expect_false(is_lognormal(normal()))
  expect_false(is_lognormal(exponential(rate = 1)))
})

# --- params ---

test_that("params.lognormal returns named vector of parameters", {
  d <- lognormal(meanlog = 1.5, sdlog = 0.7)
  p <- params(d)

  expect_named(p, c("meanlog", "sdlog"))
  expect_equal(p["meanlog"], c(meanlog = 1.5))
  expect_equal(p["sdlog"], c(sdlog = 0.7))
})

# --- mean ---

test_that("mean.lognormal returns exp(meanlog + sdlog^2/2)", {
  d <- lognormal(meanlog = 0, sdlog = 1)
  expect_equal(mean(d), exp(0 + 1 / 2))

  d2 <- lognormal(meanlog = 2, sdlog = 0.5)
  expect_equal(mean(d2), exp(2 + 0.25 / 2))

  d3 <- lognormal(meanlog = -1, sdlog = 2)
  expect_equal(mean(d3), exp(-1 + 4 / 2))
})

# --- vcov ---

test_that("vcov.lognormal returns (exp(sdlog^2)-1)*exp(2*meanlog+sdlog^2)", {
  d <- lognormal(meanlog = 0, sdlog = 1)
  expected <- (exp(1) - 1) * exp(0 + 1)
  expect_equal(vcov(d), expected)

  d2 <- lognormal(meanlog = 2, sdlog = 0.5)
  expected2 <- (exp(0.25) - 1) * exp(4 + 0.25)
  expect_equal(vcov(d2), expected2)
})

# --- dim ---

test_that("dim.lognormal returns 1 for univariate distribution", {
  d <- lognormal()
  expect_equal(dim(d), 1)
})

# --- format / print ---

test_that("format.lognormal returns correct string", {
  d <- lognormal(meanlog = 1, sdlog = 2)
  expect_equal(format(d), "Log-normal distribution (meanlog = 1, sdlog = 2)")
})

test_that("print.lognormal outputs to console", {
  d <- lognormal(meanlog = 1, sdlog = 2)
  expect_output(print(d), "Log-normal distribution")
  expect_output(print(d), "meanlog = 1")
  expect_output(print(d), "sdlog = 2")
})

test_that("print.lognormal returns object invisibly", {
  d <- lognormal()
  out <- capture.output(ret <- print(d))
  expect_identical(ret, d)
})

# --- sampler ---

test_that("sampler.lognormal returns a function that generates samples", {
  d <- lognormal(meanlog = 0, sdlog = 1)
  samp_fn <- sampler(d)

  expect_type(samp_fn, "closure")

  samples <- samp_fn(100)
  expect_length(samples, 100)

  # Samples should be positive (log-normal is on (0, Inf))
  expect_true(all(samples > 0))
})

test_that("sampler.lognormal produces samples with approximately correct mean", {
  set.seed(42)
  d <- lognormal(meanlog = 0, sdlog = 0.5)
  samp_fn <- sampler(d)
  samples <- samp_fn(10000)

  expected_mean <- exp(0 + 0.25 / 2)
  sample_mean <- sum(samples) / length(samples)
  expect_equal(sample_mean, expected_mean, tolerance = 0.05)
})

test_that("sampler.lognormal default n = 1", {
  d <- lognormal()
  samp_fn <- sampler(d)
  samples <- samp_fn()
  expect_length(samples, 1)
})

# --- density ---

test_that("density.lognormal returns correct probability density", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  pdf <- density(d)

  # Compare with dlnorm at known points
  expect_equal(pdf(0.5), dlnorm(0.5, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
  expect_equal(pdf(1), dlnorm(1, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
  expect_equal(pdf(5), dlnorm(5, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
})

test_that("density.lognormal handles log argument correctly", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  pdf <- density(d)

  expect_equal(pdf(2, log = TRUE),
               dlnorm(2, meanlog = 1, sdlog = 0.5, log = TRUE),
               tolerance = 1e-10)
})

test_that("density.lognormal returns zero for non-positive values", {
  d <- lognormal(meanlog = 0, sdlog = 1)
  pdf <- density(d)

  expect_equal(pdf(0), 0)
  expect_equal(pdf(-1), 0)
})

# --- cdf ---

test_that("cdf.lognormal returns correct cumulative distribution", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  cdf_fn <- cdf(d)

  expect_equal(cdf_fn(0), plnorm(0, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
  expect_equal(cdf_fn(1), plnorm(1, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
  expect_equal(cdf_fn(5), plnorm(5, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
})

test_that("cdf.lognormal handles log.p argument correctly", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  cdf_fn <- cdf(d)

  expect_equal(cdf_fn(2, log.p = TRUE),
               plnorm(2, meanlog = 1, sdlog = 0.5, log.p = TRUE),
               tolerance = 1e-10)
})

# --- inv_cdf ---

test_that("inv_cdf.lognormal returns correct quantiles", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  qf <- inv_cdf(d)

  expect_equal(qf(0.5), qlnorm(0.5, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
  expect_equal(qf(0.95), qlnorm(0.95, meanlog = 1, sdlog = 0.5), tolerance = 1e-10)
})

test_that("inv_cdf.lognormal round-trips with cdf", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  cdf_fn <- cdf(d)
  qf <- inv_cdf(d)

  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (p in probs) {
    expect_equal(cdf_fn(qf(p)), p, tolerance = 1e-10)
  }
})

test_that("inv_cdf.lognormal handles lower.tail and log.p", {
  d <- lognormal(meanlog = 0, sdlog = 1)
  qf <- inv_cdf(d)

  # lower.tail = FALSE means P(X > q) = p
  expect_equal(qf(0.1, lower.tail = FALSE),
               qlnorm(0.1, meanlog = 0, sdlog = 1, lower.tail = FALSE),
               tolerance = 1e-10)

  # log.p = TRUE means p is on log scale
  expect_equal(qf(log(0.5), log.p = TRUE),
               qlnorm(log(0.5), meanlog = 0, sdlog = 1, log.p = TRUE),
               tolerance = 1e-10)
})

# --- sup ---

test_that("sup.lognormal returns correct support interval", {
  d <- lognormal()
  s <- sup(d)

  expect_s3_class(s, "interval")
  expect_equal(s$infimum(), 0)
  expect_equal(s$supremum(), Inf)
  # Lower bound should be open (not including 0)
  expect_false(s$lower_closed)
})

# --- hazard ---

test_that("hazard.lognormal equals f(t)/S(t) for several test points", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  h <- hazard(d)

  test_points <- c(0.1, 0.5, 1, 2, 5, 10)
  for (t in test_points) {
    expected <- dlnorm(t, meanlog = 1, sdlog = 0.5) /
                plnorm(t, meanlog = 1, sdlog = 0.5, lower.tail = FALSE)
    expect_equal(h(t), expected, tolerance = 1e-10,
                 label = paste("hazard at t =", t))
  }
})

test_that("hazard.lognormal handles log argument correctly", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  h <- hazard(d)

  t <- 2
  expected_h <- dlnorm(t, meanlog = 1, sdlog = 0.5) /
                plnorm(t, meanlog = 1, sdlog = 0.5, lower.tail = FALSE)
  expect_equal(h(t, log = TRUE), log(expected_h), tolerance = 1e-10)
})

# --- surv ---

test_that("surv.lognormal returns 1 - cdf at test points", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  sf <- surv(d)
  cdf_fn <- cdf(d)

  test_points <- c(0, 0.1, 0.5, 1, 2, 5, 10)
  for (t in test_points) {
    expect_equal(sf(t), 1 - cdf_fn(t), tolerance = 1e-10,
                 label = paste("survival at t =", t))
  }
})

test_that("surv.lognormal handles log.p argument", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  sf <- surv(d)

  expect_equal(sf(2, log.p = TRUE),
               plnorm(2, meanlog = 1, sdlog = 0.5, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-10)
})

test_that("surv.lognormal at t=0 is 1", {
  d <- lognormal(meanlog = 1, sdlog = 0.5)
  sf <- surv(d)
  expect_equal(sf(0), 1, tolerance = 1e-10)
})

# --- Median property: median of log-normal is exp(meanlog) ---

test_that("median of lognormal is exp(meanlog)", {
  d <- lognormal(meanlog = 2, sdlog = 0.5)
  qf <- inv_cdf(d)

  expect_equal(qf(0.5), exp(2), tolerance = 1e-10)
})

# --- Vectorized operations ---

test_that("density, cdf, surv, hazard work on vectors", {
  d <- lognormal(meanlog = 0, sdlog = 1)
  pdf <- density(d)
  cdf_fn <- cdf(d)
  sf <- surv(d)
  h <- hazard(d)

  t_vec <- c(0.1, 0.5, 1, 2, 5)

  # density
  expect_equal(pdf(t_vec), dlnorm(t_vec, meanlog = 0, sdlog = 1),
               tolerance = 1e-10)

  # cdf
  expect_equal(cdf_fn(t_vec), plnorm(t_vec, meanlog = 0, sdlog = 1),
               tolerance = 1e-10)

  # survival
  expect_equal(sf(t_vec), plnorm(t_vec, meanlog = 0, sdlog = 1,
                                  lower.tail = FALSE),
               tolerance = 1e-10)

  # hazard = f/S
  expected_h <- dlnorm(t_vec, 0, 1) / plnorm(t_vec, 0, 1, lower.tail = FALSE)
  expect_equal(h(t_vec), expected_h, tolerance = 1e-10)
})
