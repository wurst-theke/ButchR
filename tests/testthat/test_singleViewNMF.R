context("Single view NMF")
library(ButchR)

ranks <- 2:4
X <- matrix(abs(rnorm(1000)), ncol = 10)
colnames(X) <- paste0("Sample", 1:10)
n_inits <- 5
nmf_exp <- run_NMF_tensor(X,
                          ranks = ranks,
                          n_initializations = n_inits,
                          extract_features = TRUE)
h <- HMatrix(nmf_exp, k=2)
w <- WMatrix(nmf_exp, k=2)


test_that("Factorization metrics", {
  expect_output(show(nmf_exp)) # default print
  # data.frame with Frob error for each init
  expect_is(FrobError(nmf_exp), "data.frame")
  expect_equal(dim(FrobError(nmf_exp)), c(n_inits, length(ranks)))
  # data.frame with OptKStats for each init
  expect_is(OptKStats(nmf_exp), "data.frame")
  expect_equal(dim(OptKStats(nmf_exp)), c(length(ranks), 10))
  # OptK
  expect_is(OptK(nmf_exp), "integer")
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

  # X without colnames and with rownames
  X <- matrix(abs(rnorm(1000)), ncol = 10)
  rownames(X) <- paste0("f", 1:100)
  x <- run_NMF_tensor(X, ranks = 2, n_initializations = 2)
  expect_is(WMatrix(x), "list")
  expect_is(WMatrix(x , k=2), "matrix")
})


test_that("Feature extraction", {
  ssf3m <- SignatureSpecificFeatures(nmf_exp, k = 3, return_all_features = TRUE)
  # feature extraction only for K>2
  expect_error(SignatureSpecificFeatures(nmf_exp, k = 2))
  expect_error(SignatureSpecificFeatures(nmf_exp, k = 5))
  expect_length(SignatureSpecificFeatures(nmf_exp), length(ranks)-1) # All K
  expect_length(SignatureSpecificFeatures(nmf_exp, k = 3), 3) # select k
  expect_is(ssf3m, "matrix")
  expect_equal(dim(ssf3m), c(nrow(X), 3)) # dim

  expect_length(SignatureFeatures(nmf_exp), length(ranks)-1) # All K
  expect_is(SignatureFeatures(nmf_exp), "list")
  expect_error(SignatureFeatures(nmf_exp , k = 2))
  expect_is(SignatureFeatures(nmf_exp , k = 3), "matrix")
  expect_error(SignatureFeatures(nmf_exp, k = 5))
})


test_that("Feature extraction after decomposition", {
  nmf_exp_nf <- run_NMF_tensor(X,
                               ranks = ranks,
                               n_initializations = n_inits,
                               extract_features = FALSE)
  expect_error(SignatureSpecificFeatures(nmf_exp_nf)) # if no feature extraction
  expect_error(SignatureFeatures(nmf_exp_nf)) # if no feature extraction
  nmf_exp_nf <- compute_SignatureFeatures(nmf_exp_nf) # compute features
  nmf_exp_nf_ssf <- SignatureSpecificFeatures(nmf_exp_nf, k = 3, return_all_features = TRUE)
  expect_is(nmf_exp_nf_ssf, "matrix")
})

test_that("Feature extraction after decomposition error K == 2", {
  nmf_exp_nf <- run_NMF_tensor(X,
                               ranks = 2,
                               n_initializations = n_inits,
                               extract_features = FALSE)
  expect_error(SignatureSpecificFeatures(nmf_exp_nf)) # if no feature extraction
  expect_error(compute_SignatureFeatures(nmf_exp_nf)) # K =2 not supported
})



test_that("W NMF normalization", {
  nmf_exp_wnorm <- normalizeW(nmf_exp)
  h_wnorm <- HMatrix(nmf_exp_wnorm, k = 2)
  w_wnorm <- WMatrix(nmf_exp_wnorm, k = 2)
  expect_is(nmf_exp_wnorm, "ButchR_NMF")
  expect_equal(sum(colSums(w_wnorm)), 2) # Culsums equal to 1
  expect_equal(dim(w_wnorm %*% h_wnorm), dim(w %*% h)) # Culsums equal to 1
})

test_that("H NMF normalization", {
  nmf_exp_hnorm <- normalizeH(nmf_exp)
  h_hnorm <- HMatrix(nmf_exp_hnorm, k = 2)
  w_hnorm <- WMatrix(nmf_exp_hnorm, k = 2)
  expect_is(nmf_exp_hnorm, "ButchR_NMF")
  expect_equal(sum(rowSums(h_hnorm)), 2) # Culsums equal to 1
  expect_equal(dim(w_hnorm %*% h_hnorm), dim(w %*% h)) # Culsums equal to 1
})

test_that("NMF regularize H", {
  nmf_exp_hreg <- regularizeH(nmf_exp)
  h_hreg <- HMatrix(nmf_exp_hreg, k = 2)
  w_hreg <- WMatrix(nmf_exp_hreg, k = 2)
  expect_is(nmf_exp_hreg, "ButchR_NMF")
  expect_equal(sum(matrixStats::rowMaxs(h_hreg)), 2) # Max fo each row 1
  expect_equal(dim(w_hreg %*% h_hreg), dim(w %*% h)) # Culsums equal to 1
})



test_that("NMF metric Plots", {
  gstat <- gg_plotKStats(nmf_exp)
  # factorization metrix
  expect_is(gstat, "ggplot")
  expect_length(unique(gstat$data$k), length(ranks)) # for all k
  # FrobError cv sumSilWidth meanSilWidth copheneticCoeff meanAmariDist
  expect_length(unique(gstat$data$Metric), 6)
  # for more metrics
  gstat <- gg_plotKStats(nmf_exp,
                         plot_vars = c("FrobError", "FrobError_min",
                                       "FrobError_mean", "FrobError_cv",
                                       "FrobError_sd", "sumSilWidth",
                                       "meanSilWidth", "copheneticCoeff",
                                       "meanAmariDist"))
  expect_length(unique(gstat$data$Metric), 9)
  # Errors
  expect_error(gg_plotKStats(X))
  expect_error(gg_plotKStats(nmf_exp, plot_vars = X))
  expect_error(gg_plotKStats(nmf_exp, plot_vars = c("a", "FrobError")))
})

test_that("NMF metric riverplot", {
  log <- capture.output({
    res <- plot(generateRiverplot(nmf_exp))
  })
  log2 <- capture.output({
    res <- plot(generateRiverplot(nmf_exp, useH = TRUE, ranks = 2:3))
  })
  # riverplot
  expect_is(log, "character") # all Ks
  expect_is(log2, "character") # selected Ks
  expect_error(generateRiverplot(nmf_exp, ranks = 4:7)) # K not in original decomp
  expect_error(generateRiverplot(nmf_exp, ranks = 2)) # expect a range
})


test_that("NMF recovery plots", {

  ggrecov <- recovery_plot(h, c(rep("A", 5), rep("B", 5)))
  # recovery plot
  expect_is(ggrecov, "ggplot")
  expect_length(unique(ggrecov$data$Annotation_level), 2) # for all annots
  expect_length(unique(ggrecov$data$SignatureID), nrow(h)) # for all k

})




