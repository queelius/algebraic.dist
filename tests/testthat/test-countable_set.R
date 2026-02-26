# Tests for the countable_set R6 class

test_that("countable_set constructor creates valid object with default lower bound", {
  cs <- countable_set$new()

  expect_s3_class(cs, "countable_set")
  expect_equal(cs$lower_bound, 0L)
})

test_that("countable_set constructor creates object with custom lower bound", {
  cs <- countable_set$new(lower = 1L)

  expect_equal(cs$lower_bound, 1L)
})

test_that("countable_set constructor coerces numeric to integer", {
  cs <- countable_set$new(lower = 5)

  expect_equal(cs$lower_bound, 5L)
  expect_type(cs$lower_bound, "integer")
})

test_that("has.countable_set accepts non-negative integers", {
  cs <- countable_set$new(lower = 0L)

  expect_true(has(cs, 0))
  expect_true(has(cs, 1))
  expect_true(has(cs, 100))
})

test_that("has.countable_set rejects negative values", {
  cs <- countable_set$new(lower = 0L)

  expect_false(has(cs, -1))
  expect_false(has(cs, -100))
})

test_that("has.countable_set rejects non-integers", {
  cs <- countable_set$new(lower = 0L)

  expect_false(has(cs, 1.5))
  expect_false(has(cs, 0.1))
  expect_false(has(cs, 2.999))
})

test_that("has.countable_set respects lower bound", {
  cs <- countable_set$new(lower = 3L)

  expect_false(has(cs, 0))
  expect_false(has(cs, 2))
  expect_true(has(cs, 3))
  expect_true(has(cs, 10))
})

test_that("has.countable_set works with vector input", {
  cs <- countable_set$new(lower = 0L)

  # all() is used internally, so a vector with all valid values returns TRUE

  expect_true(has(cs, c(0, 1, 2, 3)))
  # a vector with one invalid value returns FALSE
  expect_false(has(cs, c(0, 1, -1, 3)))
  expect_false(has(cs, c(0, 1.5, 2)))
})

test_that("infimum.countable_set returns lower_bound", {
  cs0 <- countable_set$new(lower = 0L)
  cs5 <- countable_set$new(lower = 5L)

  expect_equal(infimum(cs0), 0L)
  expect_equal(infimum(cs5), 5L)
})

test_that("supremum.countable_set returns Inf", {
  cs <- countable_set$new(lower = 0L)

  expect_equal(supremum(cs), Inf)
})

test_that("dim.countable_set returns 1", {
  cs <- countable_set$new(lower = 0L)

  expect_equal(dim(cs), 1)
})
