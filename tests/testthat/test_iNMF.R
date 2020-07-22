norm_mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
                      b = matrix(abs(rnorm(1000)), ncol = 10))
lambdas <- c(0, 0.6)
k <- 3

inmf_exp <- run_iNMF_tensor(norm_mat_list,
                            ranks = k,
                            n_initializations     = 2,
                            iterations            = 10^4,
                            convergence_threshold = 10,
                            extract_features = FALSE)


test_that("iNMF feature extraction", {
  # iNMF
  expect_output(show(inmf_exp)) # deafault print
  expect_error(SignatureSpecificFeatures(inmf_exp)) # if no feature extraction
  inmf_exp <- compute_SignatureFeatures(inmf_exp) # compute features
  inmf_exp_ssf <- SignatureSpecificFeatures(inmf_exp)
  expect_length(inmf_exp_ssf, length(norm_mat_list)) # one iNMF for each lambda
  expect_length(inmf_exp_ssf[[1]][[1]], k) # for each K
})

test_that("iNMF feature extraction error K == 2", {
  # iNMF
  inmf_exp <- run_iNMF_tensor(norm_mat_list,
                              ranks = 2,
                              n_initializations     = 2,
                              iterations            = 10^4,
                              convergence_threshold = 10,
                              extract_features = FALSE)


  expect_error(SignatureSpecificFeatures(inmf_exp)) # if no feature extraction
  expect_error(compute_SignatureFeatures(inmf_exp)) # K =2 not supported
})




#
# run_iNMF_tensor(norm_mat_list,
#                 ranks = 2:5,
#                 n_initializations     = 0,
#                 iterations            = 10^4,
#                 convergence_threshold = 40)@OptKStats
#
#
# inmf_exp <- run_iNMF_tensor(norm_mat_list,
#                                ranks = 2:5,
#                                n_initializations     = 2,
#                                iterations            = 10^4,
#                                convergence_threshold = 40)
# inmf_exp
#
# #
# #
# #
# iNMF_lambda_tuning(matrix_list           = norm_mat_list,
#                    lambdas               = 0,
#                    Output_type           = "residuals",
#                    rank                  = 3,
#                    n_initializations     = 5,
#                    iterations            = 10^4,
#                    convergence_threshold = 40,
#                    Sp                    = 0,
#                    extract_features      = FALSE)
#
#
#
# # Retrieve all the objects and extract the iNMF for the best lambda:
# inmf_tune <- iNMF_lambda_tuning(matrix_list           = norm_mat_list,
#                                 lambdas               = seq(0, 1, 0.1),
#                                 thr_cons              = 4,
#                                 Output_type           = "all",
#                                 rank                  = 9,
#                                 n_initializations     = 5,
#                                 iterations            = 10^4,
#                                 convergence_threshold = 40)
# min(inmf_tune$residuals$lambda[inmf_tune$residuals$best_lambda])
# inmf_tune$iNMF$lambda_0.2