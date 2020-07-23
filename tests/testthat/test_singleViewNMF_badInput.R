context("Bad inputs Single view NMF")
library(ButchR)

ranks <- 2:4
X <- matrix(abs(rnorm(1000)), ncol = 10)
colnames(X) <- paste0("Sample", 1:10)
n_inits <- 5



test_that("A matrix with negative values return error", {
  expect_error(run_NMF_tensor(matrix(rnorm(1000), ncol = 10),
                              ranks = 2, n_initializations = 1))
  expect_error(run_NMF_tensor(matrix(rep("a", 1000), ncol = 10),
                              ranks = 2, n_initializations = 1))
  expect_warning(run_NMF_tensor(as.data.frame(X),
                                ranks = 2, n_initializations = 1))
  expect_error(run_NMF_tensor(as.data.frame(matrix(rep("a", 1000)), ncol = 10),
                              ranks = 2, n_initializations = 1))
})

test_that("Bad ranks", {
  expect_error(run_NMF_tensor(X, ranks = "a", n_initializations = 1))
  expect_error(run_NMF_tensor(X, ranks = 20, n_initializations = 1))
  expect_error(run_NMF_tensor(X, ranks = -4:1, n_initializations = 1))
})



