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
  expect_error(run_NMF_tensor("a", ranks = 2, n_initializations = 1))
})

test_that("Bad ranks", {
  expect_error(run_NMF_tensor(X, ranks = "a", n_initializations = 1))
  expect_error(run_NMF_tensor(X, ranks = 20, n_initializations = 1)) # k greater than cols
  expect_error(run_NMF_tensor(X, ranks = -4:1, n_initializations = 1))
})



test_that("Bad method", {
  expect_error(run_NMF_tensor(X, ranks = 2, method = 1))
  expect_error(run_NMF_tensor(X, ranks = 2, method = c("nmf")))
  expect_error(run_NMF_tensor(X, ranks = 2, method = c("nmf", "nmf2")))
})



test_that("Bad n_initializations", {
  expect_error(run_NMF_tensor(X, ranks = 2, n_initializations = "a"))
  expect_error(run_NMF_tensor(X, ranks = 2, n_initializations = 1:3))
  expect_error(run_NMF_tensor(X, ranks = 2, n_initializations = 1.1))
  expect_error(run_NMF_tensor(X, ranks = 2, n_initializations = -1.1))
  expect_error(run_NMF_tensor(X, ranks = 2, n_initializations = -3))
})

test_that("Bad iterations", {
  expect_error(run_NMF_tensor(X, ranks = 2, iterations = "a"))
  expect_error(run_NMF_tensor(X, ranks = 2, iterations = 10000:10001))
  expect_error(run_NMF_tensor(X, ranks = 2, iterations = 10000.1))
  expect_error(run_NMF_tensor(X, ranks = 2, iterations = -10000.1))
  expect_error(run_NMF_tensor(X, ranks = 2, iterations = -10000))
})

test_that("Bad convergence_threshold", {
  expect_error(run_NMF_tensor(X, ranks = 2, convergence_threshold = "a"))
  expect_error(run_NMF_tensor(X, ranks = 2, convergence_threshold = 30:31))
  expect_error(run_NMF_tensor(X, ranks = 2, convergence_threshold = 30.1))
  expect_error(run_NMF_tensor(X, ranks = 2, convergence_threshold = -30.1))
  expect_error(run_NMF_tensor(X, ranks = 2, convergence_threshold = -30))
})

test_that("Bad n_neighbors", {
  expect_error(run_NMF_tensor(X, ranks = 2, n_neighbors = "a"))
  expect_error(run_NMF_tensor(X, ranks = 2, n_neighbors = 4:3))
  expect_error(run_NMF_tensor(X, ranks = 2, n_neighbors = 4.1))
  expect_error(run_NMF_tensor(X, ranks = 2, n_neighbors = -3.1))
  expect_error(run_NMF_tensor(X, ranks = 2, n_neighbors = -3))
})


test_that("Bad alpha", {
  expect_error(run_NMF_tensor(X, ranks = 2, alpha = "a"))
  expect_error(run_NMF_tensor(X, ranks = 2, alpha = 4:3))
  expect_error(run_NMF_tensor(X, ranks = 2, alpha = -3.1))
})

test_that("Bad lamb", {
  expect_error(run_NMF_tensor(X, ranks = 2, lamb = "a"))
  expect_error(run_NMF_tensor(X, ranks = 2, lamb = 4:3))
  expect_error(run_NMF_tensor(X, ranks = 2, lamb = -3.1))
})

test_that("Bad graph", {
  graph <- matrix(abs(rnorm(100)), ncol = 10)
  # Good input print
  nmf_exp <- run_NMF_tensor(X, ranks = 2, method = "GRNMF_SC", graph = graph)
  expect_output(show(nmf_exp))
  # Bad input
  expect_warning(run_NMF_tensor(X, ranks = 2, method = "NMF", graph = graph))
  expect_warning(run_NMF_tensor(X, ranks = 2, method = "GRNMF_SC", graph = as.data.frame(graph)))
  expect_error(run_NMF_tensor(X, ranks = 2, method = "GRNMF_SC", graph = "a"))
  graph <- matrix(abs(rnorm(110)), ncol = 11)
  expect_error(run_NMF_tensor(X, ranks = 2, method = "GRNMF_SC", graph = graph))
  graph <- matrix(abs(rnorm(121)), ncol = 11)
  expect_error(run_NMF_tensor(X, ranks = 2, method = "GRNMF_SC", graph = graph))
  graph <- matrix(rep("a", 100), ncol = 10)
  expect_error(run_NMF_tensor(X, ranks = 2, method = "GRNMF_SC", graph = graph))
})


test_that("Bad extract_features", {
  expect_error(run_NMF_tensor(X, ranks = 2, extract_features = "a"))
  expect_error(run_NMF_tensor(X, ranks = 2, extract_features = c(FALSE, FALSE)))
})

