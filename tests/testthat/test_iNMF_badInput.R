context("Bad inputs Integrative NMF")
library(ButchR)


mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
                 b = matrix(abs(rnorm(1000)), ncol = 10))
ranks <- 2:3
n_inits <- 2
inmf_exp <- run_iNMF_tensor(mat_list,
                            ranks = ranks,
                            n_initializations     = n_inits,
                            iterations            = 10^4,
                            convergence_threshold = 10,
                            extract_features = FALSE)



test_that("A list of matrices with negative values return error", {
  expect_error(run_iNMF_tensor(list(matrix(rnorm(1000), ncol = 10), mat_list$a),
                               ranks = 2, n_initializations = 1))
  expect_error(run_iNMF_tensor(list(matrix(rep("a", 1000), ncol = 10),
                                    mat_list$a),
                               ranks = 2, n_initializations = 1))
  expect_warning(run_iNMF_tensor(lapply(mat_list, as.data.frame),
                                ranks = 2, n_initializations = 1))
  expect_error(run_iNMF_tensor("a", ranks = 2, n_initializations = 1))
})

test_that("Non matching col names and list without names", {
  # Non matching col names
  mat_list_2 <- mat_list
  colnames(mat_list_2$a) <- paste0("sample",1:ncol(mat_list_2$a))
  expect_error(run_iNMF_tensor(mat_list_2,
                               ranks = 2, n_initializations = 1))
  # diff number of columns
  mat_list_2 <- mat_list
  mat_list_2$a <- cbind(mat_list_2$a, mat_list_2$a)
  expect_error(run_iNMF_tensor(mat_list_2,
                               ranks = 2, n_initializations = 1))
  expect_warning(run_iNMF_tensor(list(mat_list$a, mat_list$b),
                                 ranks = 2, n_initializations = 1))
})

test_that("Bad ranks", {
  expect_error(run_iNMF_tensor(mat_list, ranks = "a", n_initializations = 1))
  expect_error(run_iNMF_tensor(mat_list, ranks = 20, n_initializations = 1)) # k greater than cols
  expect_error(run_iNMF_tensor(mat_list, ranks = -4:1, n_initializations = 1))
})


test_that("Bad n_initializations", {
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, n_initializations = "a"))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, n_initializations = 1:3))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, n_initializations = 1.1))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, n_initializations = -1.1))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, n_initializations = -3))
})

test_that("Bad iterations", {
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, iterations = "a"))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, iterations = 10000:10001))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, iterations = 10000.1))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, iterations = -10000.1))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, iterations = -10000))
})

test_that("Bad convergence_threshold", {
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, convergence_threshold = "a"))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, convergence_threshold = 30:31))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, convergence_threshold = 30.1))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, convergence_threshold = -30.1))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, convergence_threshold = -30))
})

test_that("Bad Sp", {
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, Sp = "a"))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, Sp = 4:3))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, Sp = -3.1))
})

test_that("Bad lamb", {
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, lamb = "a"))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, lamb = 4:3))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, lamb = -3.1))
})

test_that("Bad extract_features", {
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, extract_features = "a"))
  expect_error(run_iNMF_tensor(mat_list, ranks = 2, extract_features = c(FALSE, FALSE)))
})

