# global reference to tensorflow (will be initialized in .onLoad)
tf <- NULL
np <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to tensorflow
  np <<- reticulate::import("numpy", delay_load = TRUE)
  tf <<- reticulate::import("tensorflow", delay_load = TRUE)


}

# Copyright © 2015-2020  The ButchR package contributors
# This file is part of the ButchR package. The ButchR package is licenced
# under GPL-3

#------------------------------------------------------------------------------#
#                         NMF GPU Wrapper - FUNCTIONS                          #
#------------------------------------------------------------------------------#
#' Single view NMF implementations
#'
#' Python functions to run several NMF implementations on tensorflow
#' using the reticulate framework, uses a non-negative matrix as input
#'
#' @param method string with method to use, "NMF", "GRNMF_SC"
#'
#' @return python function to carry out the factorization on tensorflow
#'
#' @import reticulate
#' @examples
#' \dontrun{
#' x <- ButchR:::source_NMFtensor_function("NMF")
#' x
#' }
source_NMFtensor_function <- function(method) {
  NMF_methods_list <- list(NMF      = c("nmf_tensor_lite", "NMF_tensor_py"),
                           GRNMF_SC = c("nmf_tensor_regularized_lite", "NMF_tensor_py"),
                           jNMF     = c("tensor_jNMF_lite", "jNMF_tensor_py"),
                           iNMF     = c("tensor_iNMF_lite", "iNMF_tensor_py"))
  module = NMF_methods_list[[method]]

  # Source NMF tensorflow python script
  path <- file.path(system.file(package = "ButchR"), "python/")
  tensorTheCleaver <- import_from_path("tensorTheCleaver", path = path)
  return(tensorTheCleaver[[module[1]]][[module[2]]])
}


#------------------------------------------------------------------------------#
#                      NMF tensorflow Wrapper - FUNCTION                       #
#------------------------------------------------------------------------------#

#' Single view NMF
#'
#' Computes NMF on tensorflow using the reticulate framework, uses a
#' non-negative matrix as input
#'
#' @param X Input matrix, should be a non-negative matrix
#' @param n_initializations Number of initializations to evaluate
#' @param ranks numeric vector with ranks to factorize
#' @param method method to use in the factorization, available:
#' "NMF", "GRNMF_SC"
#' @param iterations Maximum number of iterations to run for every initialization
#' @param convergence_threshold The factorization stops, if the convergence test
#'  constant for this number of iterations
#' @param n_neighbors for method "GRNMF_SC", the number of neighbors to take
#' into account when building the graph G
#' @param alpha for method "GRNMF_SC", regularization parameter alpha
#' @param lamb for method "GRNMF_SC", regularization parameter lambda
#' @param graph for method "GRNMF_SC", square matrix representing a graph
#' between columns of the input matrix, values correspond to edge weight,
#' if NULL will compute own graph.
#' @param extract_features if TRUE performs feature extraction for all f
#' actorization ranks > 2.
#'
#'
#' @return A nmfExperiment_lite object, containg the factorized matrices W and
#'  H, along with factorization metrics
#' @import reticulate
#' @export
#'
#' @examples
#' nmf_exp <- runNMFtensor_lite(matrix(1:1000, ncol = 10),
#'                              ranks = 2:5,
#'                              n_initializations     = 10,
#'                              iterations            = 10^4,
#'                              convergence_threshold = 40)
#' nmf_exp
runNMFtensor_lite <- function (X,
                               ranks,
                               method                = "NMF",
                               n_initializations     = 10,
                               iterations            = 10^4,
                               convergence_threshold = 40,
                               n_neighbors           = 4,
                               alpha                 = 0.1,
                               lamb                  = 10,
                               graph                 = NULL,
                               extract_features      = FALSE){

  #----------------------------------------------------------------------------#
  #                                    Setup                                   #
  #----------------------------------------------------------------------------#
  # Check  data
  if (min(X) < 0 ) {
    stop("\nNegative values present in input matrix\n
         only non-negative matrices supported\n")
  }
  # Convert params to integer
  nmf_params <- lapply(list(ranks                 = ranks,
                            n_initializations     = n_initializations,
                            iterations            = iterations,
                            convergence_threshold = convergence_threshold,
                            n_neighbors           = n_neighbors),
                       as.integer)
  names(nmf_params$ranks) <- paste0("k", nmf_params$ranks)

  #----------------------------------------------------------------------------#
  #                Run single view NMF on tensorflow                           #
  #----------------------------------------------------------------------------#
  # Source NMF tensorflow python function
  NMF_tensor_py <- source_NMFtensor_function(method)

  # Run NMF
  complete_eval <- lapply(nmf_params$ranks, function(k) {
    print(Sys.time())
    cat("Factorization rank: ", k, "\n")

    k_eval <- NMF_tensor_py(matrix            = X,
                            rank              = k,
                            n_initializations = nmf_params$n_initializations,
                            iterations        = nmf_params$iterations,
                            stop_threshold    = nmf_params$convergence_threshold,
                            n_neighbors       = nmf_params$n_neighbors,
                            alpha             = alpha,
                            lamb              = lamb,
                            graph             = graph)
    names(k_eval) <- c("W", "H", "iterations", "Frob_error", "W_eval")
    k_eval$iterations <- unlist(k_eval$iterations)
    k_eval$Frob_error <- unlist(k_eval$Frob_error)

    # Optimal K stats
    k_eval$OptKStats <- try(compute_OptKStats_NMF(k_eval, k), silent = FALSE)
    k_eval$W_eval <- NULL

    print(paste("NMF converged after ", paste(k_eval$iterations, collapse = ","), "iterations"))
    return(k_eval)
  })

  #----------------------------------------------------------------------------#
  #                        Build NMF object slots                              #
  #----------------------------------------------------------------------------#
  # input matrix info
  input_matrix <- list(hash       = digest::digest(X),
                       dim        = dim(X),
                       colnames   = colnames(X),
                       rownames   = rownames(X),
                       run_params = list(method            = method,
                                         n_initializations = nmf_params$n_initializations,
                                         iterations        = nmf_params$iterations,
                                         stop_threshold    = nmf_params$convergence_threshold,
                                         n_neighbors       = nmf_params$n_neighbors,
                                         alpha             = alpha,
                                         lamb              = lamb,
                                         extract_features  = extract_features))

  # Frob. error data frame
  frob_errors <- as.data.frame(do.call(cbind, lapply(complete_eval, "[[" , "Frob_error")))


  # Optimal K stats
  OptKStats <- lapply(complete_eval, "[[" , "OptKStats")
  if (!any(sapply(OptKStats, inherits, "try-error"))) {
    OptKStats <- as.data.frame(dplyr::bind_rows(OptKStats))

    # Optimal K
    indCopheneticCoeff <- which(local.maxima(OptKStats$copheneticCoeff)) # Max Cophenetic Coeff
    indMeanAmariDist   <- which(local.minima(OptKStats$meanAmariDist))   # Min Amari Dist
    OptK <- OptKStats$k[intersect(indCopheneticCoeff, indMeanAmariDist)]
    if (length(OptK) == 0) {
      #warning("No optimal K could be determined from the Optimal K stat\n")
      cat("No optimal K could be determined from the Optimal K stat\n")
    }

  } else {
    OptKStats <- data.frame()
    OptK <- integer()
    cat("Error found while computing factorization stats\nSkipping Optimal K\n")
  }

  #----------------------------------------------------------------------------#
  #                  Compute signatures specific features                      #
  #----------------------------------------------------------------------------#
  if (extract_features) {
    SignFeatures <- lapply(lapply(complete_eval, "[[" , "W"), function(W){
      if (ncol(W) == 2) {
        return(NULL)
      } else {
        rownames(W) <- rownames(X)
        return(WcomputeFeatureStats(W))
      }
    })
    SignFeatures <- data.frame(do.call(cbind, SignFeatures),
                               stringsAsFactors = FALSE)
  } else {
    SignFeatures <- data.frame()
  }


  #----------------------------------------------------------------------------#
  #                       Return nmfExperiment_lite object                     #
  #----------------------------------------------------------------------------#
  nmfExperiment_lite(input_matrix = input_matrix,
                     WMatrix      = lapply(complete_eval, "[[" , "W"),
                     HMatrix      = lapply(complete_eval, "[[" , "H"),
                     FrobError    = frob_errors,
                     OptKStats    = OptKStats,
                     OptK         = OptK,
                     SignFeatures = SignFeatures)
}

#------------------------------------------------------------------------------#
#             Criteria for optimal factorization rank - FUNCTIONS              #
#------------------------------------------------------------------------------#

#' Computes factorization optimal K stats
#'
#' @param k_eval internal object return after computing NMF with tensorflow
#' @param k factorization rank
#'
compute_OptKStats_NMF <- function(k_eval, k) {
  #----------------------------------------------------------------------------#
  #                            Frobenius error stats                           #
  #----------------------------------------------------------------------------#
  frob_errors <- k_eval[["Frob_error"]]
  min_frobError  <- min(frob_errors, na.rm = TRUE)
  sd_frobError   <- stats::sd(frob_errors, na.rm = TRUE)
  mean_frobError <- mean(frob_errors, na.rm = TRUE)
  cv_frobError   <- sd_frobError / mean_frobError

  #----------------------------------------------------------------------------#
  #     compute Silhouette Width, Cophenetic Coeff and Amari Distances         #
  #----------------------------------------------------------------------------#
  WMatrix_list <- k_eval[["W_eval"]]
  B <- length(WMatrix_list)
  concat_matrix <- do.call(cbind, WMatrix_list)
  dist_matrix   <- cosineDissMat(as.matrix(concat_matrix))

  # compute Silhouette Width
  if (length(WMatrix_list) > 1) {
    #------------------------------------------------------------------------#
    #                         compute Silhouette Width                       #
    #------------------------------------------------------------------------#
    my_pam   <- cluster::pam(dist_matrix, k = ncol(WMatrix_list[[1]]),  diss = TRUE)
    sumSilWidth  <- sum(my_pam$silinfo$widths[, "sil_width"])
    meanSilWidth <- mean(my_pam$silinfo$widths[, "sil_width"])
    #------------------------------------------------------------------------#
    #                         compute Cophenetic Coeff                       #
    #------------------------------------------------------------------------#
    # compute Cophenetic Coeff
    my_hclust <- stats::hclust(stats::as.dist(dist_matrix))
    dist_cophenetic <- as.matrix(stats::cophenetic(my_hclust))
    # take distance matrices without diagonal elements
    diag(dist_matrix) <- NA
    dist_matrix <- dist_matrix[which(!is.na(dist_matrix))]
    diag(dist_cophenetic) <- NA
    dist_cophenetic <- dist_cophenetic[which(!is.na(dist_cophenetic))]
    copheneticCoeff = unlist(stats::cor(cbind(dist_cophenetic, dist_matrix))[1, 2])

    #------------------------------------------------------------------------#
    #                         compute Amari Distances                        #
    #------------------------------------------------------------------------#
    distances_list <- unlist(lapply(1:(B - 1), function(b) {
      distances <- lapply((b + 1):B, function(b.hat) {
        amariDistance(WMatrix_list[[b]], WMatrix_list[[b.hat]])
      })
    }))
    meanAmariDist <- unlist(mean(distances_list))

  } else {
    sumSilWidth     <- NA
    meanSilWidth    <- NA
    copheneticCoeff <- NA
    meanAmariDist   <- NA
  }

  #----------------------------------------------------------------------------#
  #                         Return optimal K stats                             #
  #----------------------------------------------------------------------------#
  data.frame(rank_id = paste0("k", k),
             k       = k,
             min     = min_frobError,
             mean    = mean_frobError,
             sd      = sd_frobError,
             cv      = cv_frobError,
             sumSilWidth     = sumSilWidth,
             meanSilWidth    = meanSilWidth,
             copheneticCoeff = copheneticCoeff,
             meanAmariDist   = meanAmariDist,
             stringsAsFactors = FALSE)
}



#------------------------------------------------------------------------------#
#                       Computes Signature specific features                   #
#------------------------------------------------------------------------------#
#' Computes Signature specific features
#'
#' @param W W matrix with more than 2 signatures
#' @examples
#' \dontrun{
#' WcomputeFeatureStats(W)
#' }
WcomputeFeatureStats <- function(W) {
  #----------------------------------------------------------------------------#
  #    find features with no contribution and assign names if not present      #
  #----------------------------------------------------------------------------#
  k = ncol(W)
  if (is.null(rownames(W))) {
    rownames(W) <- paste0("f", 1:nrow(W))
  }
  idx <- rowSums(W) == 0
  # Features that contribute towards one or more signatures
  Wf <- W[!idx,]
  # Non-signature features - features with 0 contribution to all signatures
  nsf <- stats::setNames(rep(paste(rep("0", k), collapse = ""),
                             sum(idx)), rownames(W)[idx])

  #----------------------------------------------------------------------------#
  #      Run k means over all rows and assign features to the clusters         #
  #----------------------------------------------------------------------------#
  ssf <- apply(Wf, 1, function(x){
    x <- sigmoidTransform(x)
    k <- stats::kmeans(x, 2)
    max_idx <- which.max(k$centers)
    #paste(if_else(k$cluster == max_idx, "1", "0"), collapse = "")
    paste(ifelse(k$cluster == max_idx, "1", "0"), collapse = "")
  })

  # Combine with features with no contribution and sort
  sig_features <- c(ssf, nsf)
  sig_features <- sig_features[match(rownames(W), names(sig_features))]
  return(sig_features)
}





#' Create distance matrix with cosine similarity with matrix operations
#'
#' @param in.matrix concat_matrix a concatenated W matrix across multiple initializations
#' @param in.dimension factorization rank
#'
#' @examples
#' \dontrun{
#' # concat_matrix a concatenated W matrix across multiple initializations
#' cosineDissMat(as.matrix(concat_matrix))
#' }
cosineDissMat <- function(in.matrix, in.dimension=2){
  if(in.dimension == 1) in.matrix <- t(in.matrix)
  squaredVectorSum <- apply(in.matrix, 2, function(m) { sqrt(sum(m * m)) })
  squaredVectorProduct <- squaredVectorSum %*% t(squaredVectorSum)
  squaredInputSum <- t(in.matrix) %*% in.matrix
  # sum(a*b) for any a,b in M
  diss.matrix <- 1 - squaredInputSum / squaredVectorProduct
  # CosineDistance = 1 - CosineSimilarity
  return(round(diss.matrix, digits = 14))
}

#' Compute amari type distance between two matrices
#'
#' @param matrix.A,matrix.B of the same dimensionality
#'
#' @return The amari type distance of matrix.A & matrix.B according
#'        to [Wu et. al, PNAS 2016]
#'
#' @references \url{http://www.pnas.org/content/113/16/4290.long}
#' @examples
#' \dontrun{
#' amariDistance(WMatrix_list_a, WMatrix_list_b)
#' }
amariDistance <- function(matrix.A, matrix.B) {
  K <- dim(matrix.A)[2]
  C <- stats::cor(matrix.A, matrix.B)
  return(1 - (sum(apply(C, 1, max)) + sum(apply(C, 2, max))) / (2 * K))
}


local.minima <- function(x)
  ifelse(dplyr::lag(x) >= x & dplyr::lead(x) >= x, TRUE, FALSE)
local.maxima <- function(x)
  ifelse(dplyr::lag(x) <= x & dplyr::lead(x) <= x, TRUE, FALSE)

