norm_mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
                      b = matrix(abs(rnorm(1000)), ncol = 10))
lambdas <- c(0, 0.6)
k <- 3
atune_outtypes <- c("plot", "residuals", "iNMF", "all_iNMF", "all")
inmf_atune <- lapply(setNames(atune_outtypes, atune_outtypes), function(atune_outtype){
  iNMF_lambda_tuning(matrix_list           = norm_mat_list,
                     lambdas               = lambdas,
                     Output_type           = atune_outtype,
                     rank                  = k,
                     n_initializations     = 2,
                     iterations            = 10^4,
                     convergence_threshold = 10,
                     Sp                    = 0,
                     show_plot             = FALSE,
                     extract_features      = TRUE)
})


inmf_exp <- inmf_atune$iNMF
inmf_exp_ssf <- SignatureSpecificFeatures(inmf_exp)

test_that("iNMF alpha tuning", {
  # iNMF
  expect_is(inmf_exp, "integrative_NMF")
  expect_is(inmf_exp_ssf, "list")
  expect_length(inmf_exp_ssf, length(norm_mat_list)) # one iNMF for each lambda
  expect_length(inmf_exp_ssf[[1]][[1]], k) # for each K
  # Plot
  expect_length(unique(inmf_atune$plot$data$lambda), length(lambdas))
  expect_is(inmf_atune$plot, "ggplot")
  # Residuals
  expect_equal(dim(inmf_atune$residuals), c(length(lambdas), 7)) # dim
  expect_is(inmf_atune$residuals, "data.frame")
  # all_iNMF
  expect_length(inmf_atune$all_iNMF, length(lambdas)) # one iNMF for each lambda
  expect_is(inmf_atune$all_iNMF, "list")
  expect_is(inmf_atune$all_iNMF[[1]], "integrative_NMF")
  # ALL: list with iNMF jNMF NMF plot and residuals
  expect_length(inmf_atune$all, 5)
  expect_is(inmf_atune$all, "list")
})
