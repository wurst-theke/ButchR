context("Integrative NMF")
library(ButchR)


norm_mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
                      b = matrix(abs(rnorm(1000)), ncol = 10))
ranks <- 2:3
n_inits <- 2
inmf_exp <- run_iNMF_tensor(norm_mat_list,
                            ranks = ranks,
                            n_initializations     = n_inits,
                            iterations            = 10^4,
                            convergence_threshold = 10,
                            extract_features = FALSE)


test_that("iNMF print", {
  expect_is(inmf_exp, "ButchR_integrativeNMF")
  expect_output(show(inmf_exp)) # default print
  inmf_exp@OptK <- 3 # fake optK to test
  expect_output(show(inmf_exp)) # default print
  # FrobError(inmf_exp)
  expect_is(FrobError(inmf_exp), "data.frame")
  expect_equal(dim(FrobError(inmf_exp)), c(n_inits, length(ranks)))
})


test_that("H matrix", {
  expect_is(HMatrix(inmf_exp), "list")
  expect_is(HMatrix(inmf_exp, type = "shared"), "list")
  expect_length(HMatrix(inmf_exp, type = "shared"), length(ranks)) # for each K
  expect_length(HMatrix(inmf_exp, type = "viewspec"), length(ranks)) # for each K
  expect_length(HMatrix(inmf_exp, type = "total"), length(ranks)) # for each K
  # extract H matrices only for selected rank
  expect_is(HMatrix(inmf_exp, k = 2, type = "shared"), "matrix")
  expect_is(HMatrix(inmf_exp, k = 2, type = "viewspec"), "list")
  expect_is(HMatrix(inmf_exp, k = 2, type = "total"), "list")
  expect_is(HMatrix(inmf_exp, k = 2, type = "all"), "list")
  # extract H matrices only for selected view and rank
  expect_is(HMatrix(inmf_exp, k = 2, view_id = "a", type = "viewspec"), "matrix")
  expect_is(HMatrix(inmf_exp, k = 2, view_id = "a", type = "total"), "matrix")
  # errors
  expect_error(HMatrix(inmf_exp, k = 2, view_id = "c"))
  expect_error(HMatrix(inmf_exp, k = 2, type = "any"))
  expect_error(HMatrix(inmf_exp, k = 7))
})



test_that("W matrix", {
  expect_is(WMatrix(inmf_exp), "list")
  # extract W matrices only for selected rank
  expect_is(WMatrix(inmf_exp, k = 2), "list")
  expect_is(WMatrix(inmf_exp, k = 2)[[1]], "matrix")
  # extract W matrices only for selected view and rank
  expect_is(WMatrix(inmf_exp, k = 2, view_id = "a"), "matrix")
  # errors
  expect_error(WMatrix(inmf_exp, k = 2, view_id = "c"))
  expect_error(WMatrix(inmf_exp, k = 7))
})


test_that("iNMF feature extraction", {
  expect_error(SignatureSpecificFeatures(inmf_exp)) # if no feature extraction
  expect_error(SignatureFeatures(inmf_exp)) # if no feature extraction
  inmf_exp <- compute_SignatureFeatures(inmf_exp) # compute features
  inmf_exp_ssf <- SignatureSpecificFeatures(inmf_exp)
  expect_length(inmf_exp_ssf, length(norm_mat_list))
  expect_length(inmf_exp_ssf[[1]][[1]], 3) # for each K
  expect_is(SignatureSpecificFeatures(inmf_exp, view_id = "a", k = 3,
                                      return_all_features = TRUE)[[1]],
            "matrix") # specific view, all features
  # errors
  expect_error(SignatureSpecificFeatures(inmf_exp, view_id = "c")) # no view
  expect_error(SignatureSpecificFeatures(inmf_exp, k = 2))
  expect_error(SignatureSpecificFeatures(inmf_exp, k = 7))

  # SignatureFeatures
  expect_is(SignatureFeatures(inmf_exp), "list")
  expect_is(SignatureFeatures(inmf_exp , k = 3), "list")
  expect_is(SignatureFeatures(inmf_exp , k = 3, view_id = "a"), "list")
  # errors
  expect_error(SignatureFeatures(inmf_exp, view_id = "c")) # no view
  expect_error(SignatureFeatures(inmf_exp, k = 2))
  expect_error(SignatureFeatures(inmf_exp, k = 7))

})

test_that("iNMF feature extraction error K == 2", {
  inmf_exp <- run_iNMF_tensor(norm_mat_list,
                              ranks = 2,
                              n_initializations     = 2,
                              iterations            = 10^4,
                              convergence_threshold = 10,
                              extract_features = FALSE)
  expect_error(SignatureSpecificFeatures(inmf_exp)) # if no feature extraction
  expect_error(compute_SignatureFeatures(inmf_exp)) # K =2 not supported
})

