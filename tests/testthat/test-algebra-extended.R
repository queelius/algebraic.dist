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

# ---- Uniform location shift rules ----

test_that("uniform + scalar simplifies", {
  d <- uniform_dist(2, 5) + 3
  expect_true(is_uniform_dist(d))
  expect_equal(d$min, 5)
  expect_equal(d$max, 8)
})

test_that("scalar + uniform simplifies", {
  d <- 3 + uniform_dist(2, 5)
  expect_true(is_uniform_dist(d))
  expect_equal(d$min, 5)
  expect_equal(d$max, 8)
})

test_that("uniform - scalar simplifies", {
  d <- uniform_dist(2, 5) - 1
  expect_true(is_uniform_dist(d))
  expect_equal(d$min, 1)
  expect_equal(d$max, 4)
})

# ---- Uniform scale rules ----

test_that("positive scalar * uniform simplifies", {
  d <- 2 * uniform_dist(1, 4)
  expect_true(is_uniform_dist(d))
  expect_equal(d$min, 2)
  expect_equal(d$max, 8)
})

test_that("negative scalar * uniform flips bounds", {
  d <- (-2) * uniform_dist(1, 4)
  expect_true(is_uniform_dist(d))
  expect_equal(d$min, -8)
  expect_equal(d$max, -2)
})

test_that("uniform / scalar simplifies via multiplication", {
  d <- uniform_dist(2, 6) / 2
  expect_true(is_uniform_dist(d))
  expect_equal(d$min, 1)
  expect_equal(d$max, 3)
})

# ---- Unary negation of uniform ----

test_that("-uniform flips bounds", {
  d <- -uniform_dist(1, 4)
  expect_true(is_uniform_dist(d))
  expect_equal(d$min, -4)
  expect_equal(d$max, -1)
})

# ---- Weibull scale rule ----

test_that("positive scalar * weibull simplifies", {
  d <- 3 * weibull_dist(shape = 2, scale = 5)
  expect_true(is_weibull_dist(d))
  expect_equal(d$shape, 2)
  expect_equal(d$scale, 15)
})

test_that("weibull / scalar simplifies", {
  d <- weibull_dist(shape = 2, scale = 6) / 3
  expect_true(is_weibull_dist(d))
  expect_equal(d$shape, 2)
  expect_equal(d$scale, 2)
})

test_that("negative scalar * weibull stays edist", {
  d <- (-2) * weibull_dist(shape = 2, scale = 5)
  expect_true(is_edist(d))
})

# ---- ChiSq scale rule ----

test_that("positive scalar * chi_squared -> gamma", {
  d <- 3 * chi_squared(4)
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 2)    # df/2 = 4/2
  expect_equal(d$rate, 1/6)   # 1/(2*3)
})

test_that("chi_squared / scalar -> gamma", {
  d <- chi_squared(6) / 2
  expect_true(is_gamma_dist(d))
  expect_equal(d$shape, 3)    # 6/2
  expect_equal(d$rate, 1)     # 1/(2*0.5)
})

# ---- LogNormal scale rule ----

test_that("positive scalar * lognormal simplifies", {
  d <- 2 * lognormal(3, 1)
  expect_true(is_lognormal(d))
  expect_equal(d$meanlog, 3 + log(2))
  expect_equal(d$sdlog, 1)
})

test_that("lognormal / scalar simplifies", {
  d <- lognormal(3, 1) / 2
  expect_true(is_lognormal(d))
  expect_equal(d$meanlog, 3 + log(0.5))
  expect_equal(d$sdlog, 1)
})

test_that("negative scalar * lognormal stays edist", {
  d <- (-2) * lognormal(3, 1)
  expect_true(is_edist(d))
})

# ---- LogNormal * LogNormal rule ----

test_that("lognormal * lognormal simplifies", {
  d <- lognormal(1, 2) * lognormal(3, 4)
  expect_true(is_lognormal(d))
  expect_equal(d$meanlog, 4)
  expect_equal(d$sdlog, sqrt(4 + 16))
})

test_that("prod of lognormals chains correctly", {
  d <- prod(lognormal(1, 1), lognormal(2, 1), lognormal(3, 1))
  expect_true(is_lognormal(d))
  expect_equal(d$meanlog, 6)
  expect_equal(d$sdlog, sqrt(3))
})

# ---- Division operator ----

test_that("scalar / dist creates edist", {
  d <- 1 / normal(1, 1)
  expect_true(is_edist(d))
})

test_that("dist / dist creates edist", {
  d <- normal(1, 1) / normal(2, 1)
  expect_true(is_edist(d))
})

# ---- MC validation of new rules ----

test_that("uniform shift: MC moments match closed form", {
  d_unsimplified <- edist(quote(x + 3), list(x = uniform_dist(0, 1)))
  d_simplified <- uniform_dist(3, 4)

  set.seed(42)
  s1 <- sampler(d_unsimplified)(50000)
  expect_equal(mean(s1), mean(d_simplified), tolerance = 0.05)
  expect_equal(var(s1), vcov(d_simplified), tolerance = 0.05)
})

test_that("weibull scale: MC moments match closed form", {
  d_simplified <- 3 * weibull_dist(shape = 2, scale = 5)

  set.seed(42)
  # Sample from unsimplified version
  raw <- 3 * rweibull(50000, shape = 2, scale = 5)
  expect_equal(mean(raw), mean(d_simplified), tolerance = 0.2)
})

test_that("lognormal product: MC moments match closed form", {
  d_simplified <- lognormal(1, 2) * lognormal(3, 4)

  set.seed(42)
  raw <- rlnorm(50000, 1, 2) * rlnorm(50000, 3, 4)
  expect_equal(mean(log(raw)), d_simplified$meanlog, tolerance = 0.1)
})
