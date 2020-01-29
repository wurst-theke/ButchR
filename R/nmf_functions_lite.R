# global reference to tensorflow (will be initialized in .onLoad)
tf <- NULL
np <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to tensorflow
  np <<- reticulate::import("numpy", delay_load = TRUE)
  tf <<- reticulate::import("tensorflow", delay_load = TRUE)


}

# Copyright © 2015-2020  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#==============================================================================#
#                         NMF GPU Wrapper - FUNCTIONS                          #
#==============================================================================#
#' Single view NMF implementations
#'
#' Python functions to run several NMF implementations on tensorflow
#' using the reticulate framework, uses a non-negative matrix as input
#'
#' @param method string with method to use, "NMF", "GRNMF_SC"
#'
#' @return python function to carry out the factorization on tensorflow
#'
#'
#' @examples
#' x <- source_NMFtensor_function("NMF")
#' x
#' # Bratwurst:::source_NMFtensor_function("NMF")
source_NMFtensor_function <- function(method) {
  NMF_methods_list <- list(NMF      = c("nmf_tensor_lite", "NMF_tensor_py"),
                           GRNMF_SC = c("nmf_tensor_regularized_lite", "NMF_tensor_py"),
                           jNMF     = c("tensor_jNMF_lite", "jNMF_tensor_py"))
  module = NMF_methods_list[[method]]

  # Source NMF tensorflow python script
  #path <- file.path(system.file(package = "Bratwurst"), "python/tensorBratwurst")
  path <- file.path(system.file(package = "Bratwurst"), "python/")
  tensorBratwurst <- import_from_path("tensorBratwurst", path = path)
  #return(tensorBratwurst[module[1]][module[2]])
  return(tensorBratwurst[[module[1]]][[module[2]]])
}


#==============================================================================#
#                      NMF tensorflow Wrapper - FUNCTION                       #
#==============================================================================#

#' Single view NMF
#'
#' Computes NMF on tensorflow using the reticulate framework, uses a non-negative matrix as input
#'
#' @param X Input matrix, should be a non-negative matrix
#' @param n_initializations Number of initializations to evaluate
#' @param ranks numeric vector with ranks to factorize
#' @param method method to use in the factorization, available:
#' "NMF", "GRNMF_SC"
#' @param iterations Maximum number of iterations to run for every initialization
#' @param convergence_threshold The factorization stops, if the convergence test is constant for this number of iterations
#' @param n_neighbors for method "GRNMF_SC", the number of neighbors to take into account when building the graph G
#' @param alpha for method "GRNMF_SC", regularization parameter alpha
#' @param lamb for method "GRNMF_SC", regularization parameter alpha
#'
#'
#' @return A nmfExperiment_lite object, containg the factorized matrices W and H, along with factorization metrics
#'
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
                               method = "NMF",
                               n_initializations     = 10,
                               iterations            = 10^4,
                               convergence_threshold = 40,
                               n_neighbors = 4,
                               alpha = 0.1,
                               lamb  = 10 ){

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
  #source_python(file.path(system.file(package = "Bratwurst"), "python/nmf_tensor_lite.py"))

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
                            lamb              = lamb)
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
  input_matrix <- list(hash     = digest::digest(X),
                       dim      = dim(X),
                       colnames = colnames(X),
                       rownames = rownames(X))

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
  #                       Return nmfExperiment_lite object                     #
  #----------------------------------------------------------------------------#
  nmfExperiment_lite(input_matrix = input_matrix,
                     WMatrix      = lapply(complete_eval, "[[" , "W"),
                     HMatrix      = lapply(complete_eval, "[[" , "H"),
                     FrobError    = frob_errors,
                     OptKStats    = OptKStats,
                     OptK         = OptK)
                     #FeatureStats              = "data.frame",
                     #SignatureSpecificFeatures = "list" )
                     #Factorization_ranks       = Factorization_ranks)
}

#==============================================================================#
#             Criteria for optimal factorization rank - FUNCTIONS              #
#==============================================================================#

#' Computes factorization optimal K stats
#'
#' @param complete_eval internal object return after computing NMF with tensorflow
#'
#' @examples
compute_OptKStats_NMF <- function(k_eval, k) {
  #----------------------------------------------------------------------------#
  #                            Frobenius error stats                           #
  #----------------------------------------------------------------------------#
  frob_errors <- k_eval[["Frob_error"]]
  min_frobError  <- min(frob_errors, na.rm = TRUE)
  sd_frobError   <- sd(frob_errors, na.rm = TRUE)
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
    my_hclust <- hclust(as.dist(dist_matrix))
    dist_cophenetic <- as.matrix(cophenetic(my_hclust))
    # take distance matrices without diagonal elements
    diag(dist_matrix) <- NA
    dist_matrix <- dist_matrix[which(!is.na(dist_matrix))]
    diag(dist_cophenetic) <- NA
    dist_cophenetic <- dist_cophenetic[which(!is.na(dist_cophenetic))]
    copheneticCoeff = unlist(cor(cbind(dist_cophenetic, dist_matrix))[1, 2])

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


#' #==============================================================================#
#' #                         H-MATRIX ANALYSIS FUNCTIONS                          #
#' #==============================================================================#
#'
#' #' Regularize the signatures matrix (H)
#' #'
#' #' After row regularization of the matrix H, the inverse factors are
#' #' mutiplied with the columns of W in order to keep the matrix product W*H
#' #' constant.
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return A data structure of type nmfExperiment
#' #'
#' #' @export
#' #'
#' #' @examples
#' #'  NULL
#' #'
#' regularizeH <- function(nmf.exp){
#'   # account for WMatrixList and HMatrixList
#'   all_list <- lapply(seq_along(WMatrixList(nmf.exp)), function(k_ind){
#'     k_list <-
#'       lapply(seq_along(WMatrixList(nmf.exp)[[k_ind]]), function(init_ind){
#'         tempW <- WMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         tempH <- HMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         normFactor <- rowMax(tempH)
#'         newExpo <- tempH / normFactor
#'         newSigs <- tempW * normFactor
#'         return(list(W = newSigs,
#'                     H = newExpo))
#'       })
#'     names(k_list) <- names(WMatrixList(nmf.exp)[[k_ind]])
#'     return(k_list)
#'   })
#'   names(all_list) <- names(WMatrixList(nmf.exp))
#'   thisWMatrixList <- lapply(all_list, function(current_k_list){
#'     kWMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$W)
#'     })
#'   })
#'   nmf.exp <- setWMatrixList(nmf.exp, thisWMatrixList)
#'   thisHMatrixList <- lapply(all_list, function(current_k_list){
#'     kHMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$H)
#'     })
#'   })
#'   nmf.exp <- setHMatrixList(nmf.exp, thisHMatrixList)
#'   return(nmf.exp)
#' }
#'
#'
#' #' Compute signature specific features
#' #'
#' #' @param nmf.exp A NMF experiment object
#' #' @param rowDataId The index of the rowData(nmf.exp) data.frame that should be used for
#' #'  feature extraction. In case rowData(nmf.exp)[,1] is a GRanges or a related object like
#' #'  GenomicInteractions this parameter can be ignored
#' #'
#' #' @return nmf.exp with filles SignatureSpecificFeatures container
#' #'
#' #'
#' #' @export
#' #'
#' #' @examples
#' #'
#' computeSignatureSpecificFeatures <- function(nmf.exp, rowDataId = 3){
#'   if (length(OptK(nmf.exp)) == 0){
#'     stop("You need to first define an optimal k before being able to compute
#'            signature specific features!")
#'   } else {
#'     if (nrow(FeatureStats(nmf.exp)) == 0){
#'       message("Computing feature stats...")
#'       nmf.exp <- computeFeatureStats(nmf.exp)
#'     }
#'     fstats <- FeatureStats(nmf.exp)
#'     # identify unique cluster membership strings
#'     clusterMemberships <-
#'       sapply(unique(as.character(fstats$cluster)),
#'              function (x) lengths(regmatches(x, gregexpr("1", x))))
#'     sigSpecClusters <-
#'       sort(names(clusterMemberships[which(clusterMemberships == 1)]),
#'            decreasing = TRUE)
#'
#'     if (class(rowData(nmf.exp)[, 1]) %in% c("GRanges", "GInteractions",
#'                                             "GenomicInteractions")){
#'       signatureSpecificFeatures <- lapply(sigSpecClusters, function(ssc){
#'         features <- rowData(nmf.exp)[, 1][which(fstats$cluster == ssc)]
#'         return(features)
#'       })
#'     }else{
#'       signatureSpecificFeatures <- lapply(sigSpecClusters, function(ssc){
#'         features <- rowData(nmf.exp)[, rowDataId][which(fstats$cluster == ssc)]
#'         return(features)
#'       })
#'     }
#'     names(signatureSpecificFeatures) <- sigSpecClusters
#'     nmf.exp@SignatureSpecificFeatures <- signatureSpecificFeatures
#'     return(nmf.exp)
#'   }
#' }
#'
#'
