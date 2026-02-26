# Tests for extended algebraic operations and simplification rules

# ---- Scalar multiplication ----

test_that("scalar * normal simplifies", {
  d <- 3 * normal(2, 4)
  expect_true(is_normal(d))
  expect_equal(mean(d), 6)
  expect_equal(vcov(d), 36)
})

test_that("normal * scalar simplifies", {
  d <- normal(2, 4) * 3
  expect_true(is_normal(d))
  expect_equal(mean(d), 6)
  expect_equal(vcov(d), 36)
})

test_that("scalar * gamma simplifies", {
  d <- 2 * gamma_dist(3, 4)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 3)
  expect_equal(d$rate, 2)  # 4/2
})

test_that("scalar * exponential -> gamma", {
  d <- 2 * exponential(3)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 1)
  expect_equal(d$rate, 1.5)  # 3/2
})

test_that("negative scalar * gamma stays edist", {
  d <- (-2) * gamma_dist(3, 4)
  expect_true(is_edist(d))
})

test_that("dist * dist stays edist", {
  d <- normal(0, 1) * normal(0, 1)
  expect_true(is_edist(d))
})

# ---- Location shift (dist + scalar, scalar + dist, dist - scalar) ----

test_that("normal + scalar simplifies", {
  d <- normal(2, 4) + 3
  expect_true(is_normal(d))
  expect_equal(mean(d), 5)
  expect_equal(vcov(d), 4)
})

test_that("scalar + normal simplifies", {
  d <- 3 + normal(2, 4)
  expect_true(is_normal(d))
  expect_equal(mean(d), 5)
  expect_equal(vcov(d), 4)
})

test_that("normal - scalar simplifies", {
  d <- normal(5, 4) - 3
  expect_true(is_normal(d))
  expect_equal(mean(d), 2)
  expect_equal(vcov(d), 4)
})

test_that("exp + scalar stays edist for non-normal", {
  d <- exponential(1) + 5
  expect_true(is_edist(d))
})

# ---- Gamma / Exponential addition rules ----

test_that("gamma + gamma same rate simplifies", {
  d <- gamma_dist(2, 3) + gamma_dist(4, 3)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 6)
  expect_equal(d$rate, 3)
})

test_that("gamma + gamma different rate stays edist", {
  d <- gamma_dist(2, 3) + gamma_dist(4, 5)
  expect_true(is_edist(d))
})

test_that("exp + exp same rate -> gamma(2, rate)", {
  d <- exponential(3) + exponential(3)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 2)
  expect_equal(d$rate, 3)
})

test_that("exp + exp different rate stays edist", {
  d <- exponential(3) + exponential(5)
  expect_true(is_edist(d))
})

test_that("gamma + exp same rate simplifies", {
  d <- gamma_dist(3, 2) + exponential(2)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 4)
  expect_equal(d$rate, 2)
})

test_that("exp + gamma same rate simplifies (reverse order)", {
  d <- exponential(2) + gamma_dist(3, 2)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 4)
  expect_equal(d$rate, 2)
})

# ---- Chi-squared addition ----

test_that("chi_squared + chi_squared simplifies", {
  d <- chi_squared(3) + chi_squared(5)
  expect_true(is_chi_squared(d))
  expect_equal(d$df, 8)
})

# ---- Poisson addition ----

test_that("poisson + poisson simplifies", {
  d <- poisson_dist(3) + poisson_dist(5)
  expect_true(is_poisson_dist(d))
  expect_equal(d$lambda, 8)
})

# ---- Power operator ----

test_that("N(0,1)^2 -> chi_squared(1)", {
  d <- normal(0, 1) ^ 2
  expect_true(is_chi_squared(d))
  expect_equal(d$df, 1)
})

test_that("N(0,1)^2 + N(0,1)^2 -> chi_squared(2) (chaining)", {
  d <- normal(0, 1)^2 + normal(0, 1)^2
  expect_true(is_chi_squared(d))
  expect_equal(d$df, 2)
})

test_that("non-standard normal ^2 stays edist", {
  d <- normal(1, 1) ^ 2
  expect_true(is_edist(d))
})

# ---- Math group generic ----

test_that("exp(normal) -> lognormal", {
  d <- exp(normal(2, 3))
  expect_true(is_lognormal(d))
  expect_equal(d$meanlog, 2)
  expect_equal(d$sdlog, sqrt(3))
})

test_that("log(lognormal) -> normal", {
  d <- log(lognormal(2, 3))
  expect_true(is_normal(d))
  expect_equal(mean(d), 2)
  expect_equal(vcov(d), 9)  # 3^2
})

test_that("sqrt(normal) stays edist", {
  d <- sqrt(normal(4, 1))
  expect_true(is_edist(d))
})

test_that("abs(normal) stays edist", {
  d <- abs(normal(0, 1))
  expect_true(is_edist(d))
})

# ---- Summary group generic ----

test_that("sum of normals simplifies via Summary.dist", {
  d <- sum(normal(1, 1), normal(2, 3), normal(3, 5))
  expect_true(is_normal(d))
  expect_equal(mean(d), 6)
  expect_equal(vcov(d), 9)
})

test_that("min of exponentials simplifies", {
  d <- min(exponential(1), exponential(2), exponential(3))
  expect_true(is_exponential(d))
  expect_equal(d$rate, 6)
})

test_that("min of mixed stays edist", {
  d <- min(exponential(1), normal(1, 1))
  expect_true(is_edist(d))
})

test_that("max stays edist", {
  d <- max(exponential(1), exponential(2))
  expect_true(is_edist(d))
})

test_that("prod creates edist", {
  d <- prod(normal(1, 1), normal(2, 1))
  expect_true(is_edist(d))
})

test_that("single element sum returns self", {
  n <- normal(1, 1)
  d <- sum(n)
  expect_true(is_normal(d))
  expect_equal(mean(d), 1)
})

test_that("single element min returns self", {
  e <- exponential(2)
  d <- min(e)
  expect_true(is_exponential(d))
  expect_equal(d$rate, 2)
})

# ---- Sampler n=1 bug fix ----

test_that("sampler.edist works with n=1", {
  x <- normal(0, 1)
  y <- normal(0, 1)
  ed <- edist(quote(x + y), list(x = x, y = y))
  samp <- sampler(ed)(1)
  expect_type(samp, "double")
  expect_length(samp, 1)
})

# ---- Chaining / composition ----

test_that("3 * normal(0,1)^2 works via chaining", {
  # normal(0,1)^2 = chi_squared(1)
  # 3 * chi_squared(1) = ?  -- no rule for scalar * chi_squared, stays edist
  # but this tests that the chain works without errors
  d <- 3 * normal(0, 1)^2
  # normal(0,1)^2 simplifies to chi_squared(1), then 3 * chi_squared(1) is edist
  expect_true(is_edist(d) || is_chi_squared(d) || is_gamma_dist(d))
})

test_that("sum of 4 exponentials (same rate) via Summary", {
  d <- sum(exponential(2), exponential(2), exponential(2), exponential(2))
  # Reduce with + chains: ((exp+exp)+exp)+exp
  # exp(2)+exp(2) -> gamma(2,2)
  # gamma(2,2)+exp(2) -> gamma(3,2)
  # gamma(3,2)+exp(2) -> gamma(4,2)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 4)
  expect_equal(d$rate, 2)
})

test_that("zero scalar multiplication of normal", {
  d <- 0 * normal(5, 3)
  expect_true(is_normal(d))
  expect_equal(mean(d), 0)
  expect_equal(vcov(d), 0)
})

test_that("identity scalar multiplication of normal", {
  d <- 1 * normal(5, 3)
  expect_true(is_normal(d))
  expect_equal(mean(d), 5)
  expect_equal(vcov(d), 3)
})

test_that("zero location shift of normal", {
  d <- normal(5, 3) + 0
  expect_true(is_normal(d))
  expect_equal(mean(d), 5)
  expect_equal(vcov(d), 3)
})
