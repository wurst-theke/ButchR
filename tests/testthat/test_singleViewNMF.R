context("Signle view NMF")
library(ButchR)

ranks <- 2:4
X <- matrix(1:1000, ncol = 10)
nmf_exp <- runNMFtensor_lite(X, ranks = ranks, n_initializations = 5,
                             extract_features = TRUE)
h <- HMatrix(nmf_exp, k=2)
w <- WMatrix(nmf_exp, k=2)


test_that("A matrix with negative values return error", {
  expect_error(runNMFtensor_lite(matrix(-9:990, ncol = 10),
                                 ranks = 2, n_initializations = 1))
  expect_warning(runNMFtensor_lite(as.data.frame(X),
                                   ranks = 2, n_initializations = 1))
})


test_that("Results of runNMFtensor are matrices", {
  expect_length(HMatrix(nmf_exp), length(ranks)) # default return list
  expect_length(WMatrix(nmf_exp), length(ranks)) #default return list
  expect_error(HMatrix(nmf_exp, k=5)) # out of input ranks
  expect_error(WMatrix(nmf_exp, k=5)) # out of input ranks
  expect_error(HMatrix(nmf_exp, k=2:3)) # more than one K
  expect_error(WMatrix(nmf_exp, k=2:3)) # more than one K
  expect_is(h, "matrix")
  expect_is(w, "matrix")
  expect_equal(dim(h), c(2, ncol(X))) # dim H matrix
  expect_equal(dim(w), c(nrow(X), 2)) # dim W matrix
})


test_that("Feature extraction", {
  ssf3m <- SignatureSpecificFeatures(nmf_exp, k = 3, return_all_features = TRUE)

  # feature extraction only for K>2
  expect_error(SignatureSpecificFeatures(nmf_exp, k = 2))
  expect_length(SignatureSpecificFeatures(nmf_exp), length(ranks)-1) # All K
  expect_length(SignatureSpecificFeatures(nmf_exp, k = 3), 3) # select k
  SignatureSpecificFeatures(nmf_exp, k = 3, return_all_features = TRUE)
  expect_is(ssf3m, "matrix")
  expect_equal(dim(ssf3m), c(nrow(X), 3)) # dim
})

