# Tests for the finite_set R6 class

test_that("finite_set constructor creates valid object from vector", {
  # Given: A numeric vector
  values <- c(1, 2, 3, 4, 5)

  # When: Creating a finite_set
  fs <- finite_set$new(values)

  # Then: Object has correct class and stores unique values
  expect_s3_class(fs, "finite_set")
  expect_equal(fs$values, values)
})

test_that("finite_set constructor stores unique values only", {
  # Given: A vector with duplicates
  values <- c(1, 1, 2, 2, 3)

  # When: Creating a finite_set
  fs <- finite_set$new(values)

  # Then: Only unique values are stored
  expect_equal(sort(fs$values), c(1, 2, 3))
})

test_that("finite_set works with matrix input", {
  # Given: A matrix of values
  values <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)

  # When: Creating a finite_set
  fs <- finite_set$new(values)

  # Then: Object stores the matrix
  expect_true(is.matrix(fs$values))
})

test_that("finite_set$has correctly checks membership for vectors", {
  fs <- finite_set$new(c(1, 2, 3))

  expect_true(fs$has(1))
  expect_true(fs$has(2))
  expect_true(fs$has(3))
  expect_false(fs$has(4))
  expect_false(fs$has(0))
})

test_that("finite_set$has correctly checks membership for matrices", {
  values <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  fs <- finite_set$new(values)

  # Row that exists in the set
  expect_true(fs$has(c(1, 4)))
  expect_true(fs$has(c(2, 5)))

  # Row that does not exist
  expect_false(fs$has(c(1, 5)))
  expect_false(fs$has(c(10, 20)))
})

test_that("finite_set$infimum returns minimum for vectors", {
  fs <- finite_set$new(c(3, 1, 4, 1, 5, 9, 2, 6))

  expect_equal(fs$infimum(), 1)
})

test_that("finite_set$infimum returns column-wise minimums for matrices", {
  values <- matrix(c(1, 2, 3, 10, 5, 15), nrow = 3, ncol = 2)
  fs <- finite_set$new(values)

  expect_equal(fs$infimum(), c(1, 5))
})

test_that("finite_set$supremum returns maximum for vectors", {
  fs <- finite_set$new(c(3, 1, 4, 1, 5, 9, 2, 6))

  expect_equal(fs$supremum(), 9)
})

test_that("finite_set$supremum returns column-wise maximums for matrices", {
  values <- matrix(c(1, 2, 3, 10, 5, 15), nrow = 3, ncol = 2)
  fs <- finite_set$new(values)

  expect_equal(fs$supremum(), c(3, 15))
})

test_that("finite_set$dim returns 1 for vectors", {
  fs <- finite_set$new(c(1, 2, 3))

  expect_equal(fs$dim(), 1)
})

test_that("finite_set$dim returns number of columns for matrices", {
  values <- matrix(1:12, nrow = 4, ncol = 3)
  fs <- finite_set$new(values)

  expect_equal(fs$dim(), 3)
})

test_that("has.finite_set S3 wrapper works correctly", {
  fs <- finite_set$new(c(1, 2, 3))

  expect_true(has(fs, 1))
  expect_false(has(fs, 4))
})

test_that("infimum.finite_set S3 wrapper works correctly", {
  fs <- finite_set$new(c(5, 3, 7))

  expect_equal(infimum(fs), 3)
})

test_that("supremum.finite_set S3 wrapper works correctly", {
  fs <- finite_set$new(c(5, 3, 7))

  expect_equal(supremum(fs), 7)
})

test_that("dim.finite_set S3 wrapper works correctly", {
  fs <- finite_set$new(c(1, 2, 3))

  expect_equal(dim(fs), 1)
})

test_that("finite_set handles single value", {
  fs <- finite_set$new(42)

  expect_equal(fs$values, 42)
  expect_true(fs$has(42))
  expect_false(fs$has(0))
  expect_equal(fs$infimum(), 42)
  expect_equal(fs$supremum(), 42)
  expect_equal(fs$dim(), 1)
})

test_that("finite_set handles negative values", {
  fs <- finite_set$new(c(-5, -3, 0, 2, 4))

  expect_true(fs$has(-5))
  expect_true(fs$has(0))
  expect_equal(fs$infimum(), -5)
  expect_equal(fs$supremum(), 4)
})
