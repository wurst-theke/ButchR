context("Aux functions")
library(ButchR)

rand_mat <- matrix(abs(rnorm(1000)), ncol = 10)



test_that("normalizeUpperQuartile", {
  x <- normalizeUpperQuartile(rand_mat)
  expect_is(x, "matrix")
  expect_equal(dim(x), dim(rand_mat))
})

test_that("sigmoidTransform2", {
  x <- sigmoidTransform2(rand_mat[1,])
  expect_is(x, "numeric")
  expect_length(x, ncol(rand_mat))
})

test_that("rankTransform", {
  x <- rankTransform(rand_mat)
  expect_is(x, "matrix")
  expect_equal(dim(x), dim(rand_mat))
})

test_that("orderBinary", {
  x <- orderBinary(rand_mat)
  expect_is(x, "integer")
  expect_length(x, ncol(rand_mat))
})

