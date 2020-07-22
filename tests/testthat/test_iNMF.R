norm_mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
                      b = matrix(abs(rnorm(1000)), ncol = 10))

#' # Compare the residual across multiple lambdas:
inmf_atune <- iNMF_lambda_tuning(matrix_list           = norm_mat_list,
                                 lambdas               = c(0, 0.6),
                                 Output_type           = "all",
                                 rank                  = 3,
                                 n_initializations     = 2,
                                 iterations            = 10^4,
                                 convergence_threshold = 10,
                                 Sp                    = 0,
                                 show_plot             = FALSE,
                                 extract_features      = FALSE)




test_that("iNMF alpha tuning", {
  expect_length(inmf_atune, 5) # list with iNMF jNMF NMF plot and residuals
  expect_is(inmf_atune, "list")
})






# inmf_exp <- run_iNMF_tensor(norm_mat_list,
#                                ranks = 2:5,
#                                n_initializations     = 2,
#                                iterations            = 10^4,
#                                convergence_threshold = 40)
# inmf_exp

#
#
#
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