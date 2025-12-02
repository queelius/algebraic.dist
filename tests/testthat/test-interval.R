# Tests for the interval R6 class

test_that("interval constructor creates valid object with defaults", {
  # Given: Default constructor call
  # When: Creating an interval
  i <- interval$new()

  # Then: Object represents the entire real line
  expect_s3_class(i, "interval")
  expect_equal(i$lower, -Inf)
  expect_equal(i$upper, Inf)
  expect_false(i$lower_closed)
  expect_false(i$upper_closed)
})

test_that("interval constructor creates bounded interval", {
  # Given: Specific bounds
  i <- interval$new(lower = 0, upper = 10, lower_closed = TRUE, upper_closed = FALSE)

  # Then: Bounds are set correctly
  expect_equal(i$lower, 0)
  expect_equal(i$upper, 10)
  expect_true(i$lower_closed)
  expect_false(i$upper_closed)
})

test_that("interval constructor handles vector bounds for multivariate", {
  # Given: Vector bounds (for multivariate support)
  lower <- c(0, -1, -Inf)
  upper <- c(10, 1, Inf)

  i <- interval$new(lower = lower, upper = upper)

  expect_equal(i$lower, lower)
  expect_equal(i$upper, upper)
  expect_equal(i$dim(), 3)
})

test_that("interval constructor replicates bounds to match lengths", {
  # Given: Mismatched length bounds
  i <- interval$new(lower = c(0, 0), upper = 10)

  # Then: upper should be replicated
  expect_equal(i$lower, c(0, 0))
  expect_equal(i$upper, c(10, 10))
})

test_that("interval constructor replicates closed flags", {
  i <- interval$new(lower = c(0, 0), upper = c(1, 1),
                    lower_closed = TRUE, upper_closed = FALSE)

  expect_equal(i$lower_closed, c(TRUE, TRUE))
  expect_equal(i$upper_closed, c(FALSE, FALSE))
})

test_that("interval$has correctly checks membership for closed interval", {
  # [0, 10]
  i <- interval$new(lower = 0, upper = 10, lower_closed = TRUE, upper_closed = TRUE)

  expect_true(i$has(0))
  expect_true(i$has(5))
  expect_true(i$has(10))
  expect_false(i$has(-1))
  expect_false(i$has(11))
})

test_that("interval$has correctly checks membership for open interval", {
  # (0, 10)
  i <- interval$new(lower = 0, upper = 10, lower_closed = FALSE, upper_closed = FALSE)

  expect_false(i$has(0))
  expect_true(i$has(5))
  expect_false(i$has(10))
  expect_true(i$has(0.001))
  expect_true(i$has(9.999))
})

test_that("interval$has correctly checks membership for half-open interval", {
  # [0, 10)
  i <- interval$new(lower = 0, upper = 10, lower_closed = TRUE, upper_closed = FALSE)

  expect_true(i$has(0))
  expect_true(i$has(5))
  expect_false(i$has(10))
})

test_that("interval$has handles infinite bounds correctly", {
  # (-Inf, Inf) - the real line
  i <- interval$new()

  expect_true(i$has(0))
  expect_true(i$has(-1e100))
  expect_true(i$has(1e100))
})

test_that("interval$has handles vector input for multivariate", {
  i <- interval$new(lower = c(0, 0), upper = c(1, 1),
                    lower_closed = TRUE, upper_closed = TRUE)

  expect_true(all(i$has(c(0.5, 0.5))))
  expect_true(all(i$has(c(0, 0))))
  expect_true(all(i$has(c(1, 1))))
  expect_false(all(i$has(c(0.5, 2))))
})

test_that("interval$infimum returns lower bound", {
  i <- interval$new(lower = 5, upper = 10)

  expect_equal(i$infimum(), 5)
})

test_that("interval$infimum returns vector for multivariate", {
  i <- interval$new(lower = c(1, 2, 3), upper = c(4, 5, 6))

  expect_equal(i$infimum(), c(1, 2, 3))
})

test_that("interval$supremum returns upper bound", {
  i <- interval$new(lower = 5, upper = 10)

  expect_equal(i$supremum(), 10)
})

test_that("interval$supremum returns vector for multivariate", {
  i <- interval$new(lower = c(1, 2, 3), upper = c(4, 5, 6))

  expect_equal(i$supremum(), c(4, 5, 6))
})

test_that("interval$dim returns 1 for univariate", {
  i <- interval$new(lower = 0, upper = 1)

  expect_equal(i$dim(), 1)
})

test_that("interval$dim returns correct dimension for multivariate", {
  i <- interval$new(lower = rep(0, 5), upper = rep(1, 5))

  expect_equal(i$dim(), 5)
})

test_that("interval$is_empty detects empty intervals", {
  # Empty: lower > upper
  i1 <- interval$new(lower = 10, upper = 5)
  expect_true(i1$is_empty())

  # Empty: same bound but not both closed
  i2 <- interval$new(lower = 5, upper = 5, lower_closed = TRUE, upper_closed = FALSE)
  expect_true(i2$is_empty())

  # Not empty: same bound and both closed (single point)
  i3 <- interval$new(lower = 5, upper = 5, lower_closed = TRUE, upper_closed = TRUE)
  expect_false(i3$is_empty())
})

test_that("has.interval S3 wrapper works correctly", {
  i <- interval$new(lower = 0, upper = 1, lower_closed = TRUE, upper_closed = TRUE)

  expect_true(has(i, 0.5))
  expect_false(has(i, 2))
})

test_that("infimum.interval S3 wrapper works correctly", {
  i <- interval$new(lower = 3, upper = 7)

  expect_equal(infimum(i), 3)
})

test_that("supremum.interval S3 wrapper works correctly", {
  i <- interval$new(lower = 3, upper = 7)

  expect_equal(supremum(i), 7)
})

test_that("dim.interval S3 wrapper works correctly", {
  i <- interval$new(lower = c(0, 0), upper = c(1, 1))

  expect_equal(dim(i), 2)
})

test_that("print.interval produces output without error", {
  i <- interval$new(lower = 0, upper = 1, lower_closed = TRUE, upper_closed = FALSE)

  expect_output(print(i), "\\[0, 1\\)")
})

test_that("print.interval shows correct brackets for different closedness", {
  # Open interval
  i1 <- interval$new(lower = 0, upper = 1, lower_closed = FALSE, upper_closed = FALSE)
  expect_output(print(i1), "\\(0, 1\\)")

  # Closed interval
  i2 <- interval$new(lower = 0, upper = 1, lower_closed = TRUE, upper_closed = TRUE)
  expect_output(print(i2), "\\[0, 1\\]")
})

test_that("interval handles semi-infinite intervals", {
  # [0, Inf)
  i1 <- interval$new(lower = 0, upper = Inf, lower_closed = TRUE, upper_closed = FALSE)
  expect_true(i1$has(0))
  expect_true(i1$has(1e100))
  expect_false(i1$has(-1))

  # (-Inf, 0]
  i2 <- interval$new(lower = -Inf, upper = 0, lower_closed = FALSE, upper_closed = TRUE)
  expect_true(i2$has(0))
  expect_true(i2$has(-1e100))
  expect_false(i2$has(1))
})
