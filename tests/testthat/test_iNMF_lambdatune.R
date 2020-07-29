context("Integrative NMF lambda tune")
library(ButchR)

mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
                 b = matrix(abs(rnorm(1000)), ncol = 10))
lambdas <- c(0, 0.6)
k <- 3
atune_outtypes <- c("plot", "residuals", "iNMF", "all_iNMF", "all")
inmf_atune <- lapply(setNames(atune_outtypes, atune_outtypes), function(atune_outtype){
  iNMF_lambda_tuning(matrix_list           = mat_list,
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

test_that("iNMF lambda tuning", {
  # iNMF
  expect_is(inmf_exp, "ButchR_integrativeNMF")
  expect_is(inmf_exp_ssf, "list")
  expect_length(inmf_exp_ssf, length(mat_list)) # one iNMF for each lambda
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
  expect_is(inmf_atune$all_iNMF[[1]], "ButchR_integrativeNMF")
  # ALL: list with iNMF jNMF NMF plot and residuals
  expect_length(inmf_atune$all, 5)
  expect_is(inmf_atune$all, "list")
})



#------------------------------------------------------------------------------#
#                                  Bad inputs                                  #
#------------------------------------------------------------------------------#

test_that("A list of matrices with negative values return error", {
  expect_error(
    iNMF_lambda_tuning(mat_list$a, # noly one matrix returns error
                       n_initializations = 1,
                       convergence_threshold = 1)
  )
  expect_error(
    iNMF_lambda_tuning(list(matrix(rnorm(1000), ncol = 10), mat_list$a),
                       n_initializations = 1,
                       convergence_threshold = 1)
  )
  expect_error(
    iNMF_lambda_tuning(list(matrix(rep("a", 1000), ncol = 10), mat_list$a),
                       n_initializations = 1,
                       convergence_threshold = 1)
  )
})

test_that("Non matching col names and list without names", {
  # Non matching col names
  mat_list_2 <- mat_list
  colnames(mat_list_2$a) <- paste0("sample",1:ncol(mat_list_2$a))
  expect_error(
    iNMF_lambda_tuning(mat_list_2,
                       n_initializations = 1,
                       convergence_threshold = 1)
  )
  # diff number of columns
  mat_list_2 <- mat_list
  mat_list_2$a <- cbind(mat_list_2$a, mat_list_2$a)
  expect_error(
    iNMF_lambda_tuning(mat_list_2,
                       n_initializations = 1,
                       convergence_threshold = 1)
  )
  expect_warning(
    iNMF_lambda_tuning(list(mat_list$a, mat_list$b), # no names
                       n_initializations = 1,
                       convergence_threshold = 1)
  )
})

test_that("Bad lambdas", {
  expect_error(iNMF_lambda_tuning(mat_list, lambdas = "a"))
  expect_error(iNMF_lambda_tuning(mat_list, lambdas = numeric()))
  expect_error(iNMF_lambda_tuning(mat_list, lambdas = c(-3.1, 0)))
})

test_that("Bad Output_type", {
  expect_error(iNMF_lambda_tuning(mat_list, Output_type = "a"))
  expect_error(iNMF_lambda_tuning(mat_list, Output_type = 1))
  expect_error(iNMF_lambda_tuning(mat_list, Output_type = c("residuals", "all")))
})

test_that("Bad thr_cons", {
  expect_error(iNMF_lambda_tuning(mat_list, thr_cons = "a"))
  expect_error(iNMF_lambda_tuning(mat_list, thr_cons = 4:3))
  expect_error(iNMF_lambda_tuning(mat_list, thr_cons = -3.1))
})

# test_that("Bad Sp", {
#   expect_error(iNMF_lambda_tuning(mat_list, Sp = "a"))
#   expect_error(iNMF_lambda_tuning(mat_list, Sp = 4:3))
#   expect_error(iNMF_lambda_tuning(mat_list, Sp = -3.1))
# })




# test_that("Bad ranks", {
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = "a", n_initializations = 1))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 20, n_initializations = 1)) # k greater than cols
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = -4:1, n_initializations = 1))
# })
#
#
# test_that("Bad n_initializations", {
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, n_initializations = "a"))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, n_initializations = 1:3))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, n_initializations = 1.1))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, n_initializations = -1.1))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, n_initializations = -3))
# })
#
# test_that("Bad iterations", {
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, iterations = "a"))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, iterations = 10000:10001))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, iterations = 10000.1))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, iterations = -10000.1))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, iterations = -10000))
# })
#
# test_that("Bad convergence_threshold", {
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, convergence_threshold = "a"))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, convergence_threshold = 30:31))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, convergence_threshold = 30.1))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, convergence_threshold = -30.1))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, convergence_threshold = -30))
# })
#
# test_that("Bad Sp", {
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, Sp = "a"))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, Sp = 4:3))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, Sp = -3.1))
# })
#
# test_that("Bad lamb", {
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, lamb = "a"))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, lamb = 4:3))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, lamb = -3.1))
# })
#
# test_that("Bad extract_features", {
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, extract_features = "a"))
#   expect_error(iNMF_lambda_tuning(mat_list, ranks = 2, extract_features = c(FALSE, FALSE)))
# })
