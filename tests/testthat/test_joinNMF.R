context("Join NMF")
library(ButchR)


norm_mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
                      b = matrix(abs(rnorm(1000)), ncol = 10))
ranks <- 2:3
n_inits <- 2
jnmf_exp <- run_joinNMF_tensor(norm_mat_list,
                               ranks = ranks,
                               n_initializations     = n_inits,
                               iterations            = 10^4,
                               convergence_threshold = 10,
                               extract_features = FALSE)



#gg_plotKStats(jnmf_exp)

test_that("jNMF print", {
  expect_is(jnmf_exp, "ButchR_joinNMF")
  expect_output(show(jnmf_exp)) # default print
  jnmf_exp@OptK <- 3 # fake optK to test
  expect_output(show(jnmf_exp)) # default print
  # FrobError(jnmf_exp)
  expect_is(FrobError(jnmf_exp), "data.frame")
  expect_equal(dim(FrobError(jnmf_exp)), c(n_inits, length(ranks)))
})


test_that("H matrix", {
  expect_is(HMatrix(jnmf_exp), "list")
  # extract H matrices only for selected rank
  expect_is(HMatrix(jnmf_exp, k = 2), "matrix")
  # errors
  expect_error(HMatrix(jnmf_exp, k = 7))
})



test_that("W matrix", {
  expect_is(WMatrix(jnmf_exp), "list")
  # extract W matrices only for selected rank
  expect_is(WMatrix(jnmf_exp, k = 2), "list")
  expect_is(WMatrix(jnmf_exp, k = 2)[[1]], "matrix")
  # extract W matrices only for selected view and rank
  expect_is(WMatrix(jnmf_exp, k = 2, view_id = "a"), "matrix")
  # errors
  expect_error(WMatrix(jnmf_exp, k = 2, view_id = "c"))
  expect_error(WMatrix(jnmf_exp, k = 7))
})


test_that("iNMF feature extraction", {
  expect_error(SignatureSpecificFeatures(jnmf_exp)) # if no feature extraction
  jnmf_exp <- compute_SignatureFeatures(jnmf_exp) # compute features
  jnmf_exp_ssf <- SignatureSpecificFeatures(jnmf_exp)
  expect_length(jnmf_exp_ssf, length(norm_mat_list))
  expect_length(jnmf_exp_ssf[[1]][[1]], 3) # for each K
  expect_is(SignatureSpecificFeatures(jnmf_exp, view_id = "a", k = 3,
                                      return_all_features = TRUE)[[1]],
            "matrix") # specific view, all features
  # errors
  expect_error(SignatureSpecificFeatures(jnmf_exp, view_id = "c")) # no view
  expect_error(SignatureSpecificFeatures(jnmf_exp, k = 2))
  expect_error(SignatureSpecificFeatures(jnmf_exp, k = 7))
})

test_that("iNMF feature extraction error K == 2", {
  jnmf_exp <- run_joinNMF_tensor(norm_mat_list,
                                 ranks = 2,
                                 n_initializations     = 2,
                                 iterations            = 10^4,
                                 convergence_threshold = 10,
                                 extract_features = FALSE)
  expect_error(SignatureSpecificFeatures(jnmf_exp)) # if no feature extraction
  expect_error(compute_SignatureFeatures(jnmf_exp)) # K =2 not supported
})

