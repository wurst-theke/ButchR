#' @include generics.R
NULL

# Copyright © 2015-2020  The ButchR package contributors
# This file is part of the ButchR package. The ButchR package is licenced
# under GPL-3

#' The ButchR_NMF class
#'
#' @slot input_matrix list of matrices
#' @slot WMatrix list.
#' @slot HMatrix list.
#' @slot FrobError DataFrame.
#' @slot OptKStats DataFrame.
#' @slot OptK numeric.
#' @slot SignFeatures DataFrame or list.
#'
#' @return An object of class ButchR_NMF.
#' @import methods
#' @export
#'
#' @examples
#' \dontrun{
#' ButchR_NMF(input_matrix = input_matrix,
#'            WMatrix      = lapply(complete_eval, "[[" , "W"),
#'            HMatrix      = lapply(complete_eval, "[[" , "H"),
#'            FrobError    = frob_errors,
#'            OptKStats    = OptKStats,
#'            OptK         = OptK,
#'            SignFeatures = SignFeatures)
#' }
ButchR_NMF <- setClass(
  Class = "ButchR_NMF",
  slots = list(input_matrix = "list",
               HMatrix      = "list",
               WMatrix      = "list",
               FrobError    = "data.frame",
               OptKStats    = "data.frame",
               OptK         = "numeric",
               SignFeatures = "data.frame" )

)

setMethod("show",
          "ButchR_NMF",
          function(object) {
            cat("class: ButchR_NMF \n")
            #cat("Best factorization index: ", object@best_factorization_idx, "\n")
            cat("Original matrix dimension: ", object@input_matrix$dim, "\n")
            cat("Factorization performed for ranks: ", object@OptKStats$k, "\n")
            if (length(object@OptK) == 0) {
              OptK <- "Please select manualy\n"
            } else {
              OptK <- object@OptK
            }

            cat("Optimal K based on factorization metrics: ", OptK, "\n")
            rparams <- object@input_matrix$run_params
            cat("Running parameters: \n")
            cat("method = ", rparams$method, " \n")
            cat("n_initializations = ", rparams$n_initializations, " \n")
            cat("iterations = ", rparams$iterations, " \n")
            cat("stop_threshold = ", rparams$stop_threshold, " \n")
            if (rparams$method == "GRNMF_SC") {
              cat("n_neighbors = ", rparams$n_neighbors, " \n")
              cat("alpha = ", rparams$alpha, " \n")
              cat("lamb = ", rparams$lamb, " \n")
            }
            cat("extract_features = ", rparams$extract_features, " \n")

          }
)

#------------------------------------------------------------------------------#
#                             H and W matrices                                 #
#------------------------------------------------------------------------------#
#### H-Matrix (H-Matrix with smallest frobError)
#' @rdname HMatrix-methods
#' @aliases HMatrix,ANY-method
#' @export
#'
#' @examples
#' HMatrix(nmf_exp)
#' HMatrix(nmf_exp, k = 2)
setMethod("HMatrix",
          "ButchR_NMF",
          function(x, k = NULL, ...) {
            if(is.null(k)) {
              H <- x@HMatrix
              if (!is.null(x@input_matrix$colnames)) {
                H <- lapply(H, function(hi) {
                  colnames(hi) <- x@input_matrix$colnames
                  hi
                })
              }
            } else {
              if (length(k) > 1 ) {
                stop("k:", paste0(k, collapse = ","), " expecting single value",
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","),
                     "\nOr NULL to retur a list of H matrices for all Ks")
              }
              idx <- as.character(x@OptKStats$rank_id[x@OptKStats$k == k])
              if (length(idx) == 0 ) {
                stop("No H matrix present for k = ", k,
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","))
              }

              H <- x@HMatrix[[idx]]
              if (!is.null(x@input_matrix$colnames)) {
                colnames(H) <- x@input_matrix$colnames
              }
            }
            return(H)
          }
)


# W-Matrix (W-Matrix with smallest frobError)
#' @rdname WMatrix-methods
#' @aliases WMatrix,ANY-method
#' @export
#'
#' @examples
#' # For ButchR_NMF objects:
#' WMatrix(nmf_exp)
#' WMatrix(nmf_exp, k = 2)
setMethod("WMatrix",
          "ButchR_NMF",
          function(x, k = NULL, ...) {
            if(is.null(k)) {
              W <- x@WMatrix
              if (!is.null(x@input_matrix$rownames)) {
                W <- lapply(W, function(wi) {
                  rownames(wi) <- x@input_matrix$rownames
                  wi
                })
              }
            } else {
              if (length(k) > 1 ) {
                stop("k:", paste0(k, collapse = ","), " expecting single value",
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","),
                     "\nOr NULL to retur a list of W matrices for all Ks")
              }
              idx <- as.character(x@OptKStats$rank_id[x@OptKStats$k == k])
              if (length(idx) == 0 ) {
                stop("No W matrix present for k = ", k,
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","))
              }

              W <- x@WMatrix[[idx]]
              if (!is.null(x@input_matrix$rownames)) {
                rownames(W) <- x@input_matrix$rownames
              }
            }
            return(W)
          }
)

#------------------------------------------------------------------------------#
#                              Optimal K Statistics                            #
#------------------------------------------------------------------------------#
# Return Frobenius Error from all initializations

#' @rdname FrobError-methods
#' @aliases FrobError,ANY-method
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 2)
#' FrobError(nmf_exp)
#' }
setMethod("FrobError", "ButchR_NMF", function(x, ...) x@FrobError)


#### Optimal K Statistics
#' Return optimal factorization rank (K) Statistics
#'
#' @param x an object of class ButchR_NMF.
#' @param ... additional parameters.
#'
#' @return optimal K Statistics
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 2)
#' OptKStats(nmf_exp)
#' }
setMethod("OptKStats", "ButchR_NMF", function(x, ...) x@OptKStats)

#### Optimal K
#' Return optimal K
#'
#' @param x an object of class ButchR_NMF.
#' @param ... additional parameters.
#'
#' @return numeric - optimal K
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 10)
#' OptK(nmf_exp)
#' }
setMethod("OptK", "ButchR_NMF", function(x, ...) x@OptK)

#------------------------------------------------------------------------------#
#                               NMF Normalization                              #
#------------------------------------------------------------------------------#
#' Normalize matrix by columns
#'
#' @param in_matrix  matrix to normalize.
#' @param in_dimension normalize by rows == 1 or columns == 1 .
#' @return matrix normalized by columns, colsum == 1.
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 10)
#' normalize_matrix_per_dim(WMatrix(nmf_exp, k=2), 2)
#' }
normalize_matrix_per_dim <- function (in_matrix, in_dimension) {
  if (in_dimension == 1) {
    choice_ind <- which(rowSums(in_matrix) > 0)
    out_df <- matrix(0, nrow = nrow(in_matrix), ncol = ncol(in_matrix))
    out_df[choice_ind, ] <- in_matrix[choice_ind, ]/rowSums(in_matrix)[choice_ind]
  }
  else if (in_dimension == 2) {
    t_df <- t(in_matrix)
    # Initialize a 0 Matrix in case one factor is only 0s
    choice_ind <- which(rowSums(t_df) > 0)
    temp_df <- matrix(0, nrow = nrow(t_df), ncol = ncol(t_df))
    # Normalize
    temp_df[choice_ind, ] <- t_df[choice_ind, ]/rowSums(t_df)[choice_ind]
    out_df <- t(temp_df)
  }
  # assign columns and rownames
  colnames(out_df) <- colnames(in_matrix)
  rownames(out_df) <- rownames(in_matrix)
  return(out_df)
}


#' @rdname normalizeW-methods
#' @aliases normalizeW,ANY,ANY-method
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 10)
#' normalizeW(nmf_exp)
#' }
setMethod("normalizeW",
          "ButchR_NMF",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            #nmf_exp@WMatrix
            all_list <- lapply(seq_along(nmf_exp@WMatrix), function(k_ind){
              tempW <- nmf_exp@WMatrix[[k_ind]]
              tempH <- nmf_exp@HMatrix[[k_ind]]
              normFactor <- colSums(tempW)
              newSigs <- normalize_matrix_per_dim(tempW, in_dimension = 2)
              #newSigs <- as.matrix(YAPSA::normalize_df_per_dim(tempW, 2))
              newExpo <- tempH * normFactor
              return(list(W = newSigs,
                          H = newExpo))
            })
            names(all_list) <- names(nmf_exp@WMatrix)

            # Set new matrices
            nmf_exp@WMatrix <- lapply(all_list, "[[" , "W")
            nmf_exp@HMatrix <- lapply(all_list, "[[" , "H")
            return(nmf_exp)
          }
)

#' @rdname normalizeH-methods
#' @aliases normalizeH,ANY,ANY-method
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 10)
#' normalizeH(nmf_exp)
#' }
setMethod("normalizeH",
          "ButchR_NMF",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            #nmf_exp@WMatrix
            all_list <- lapply(seq_along(nmf_exp@WMatrix), function(k_ind){
              tempW <- nmf_exp@WMatrix[[k_ind]]
              tempH <- nmf_exp@HMatrix[[k_ind]]
              normFactor <- rowSums(tempH)
              newExpo <- normalize_matrix_per_dim(tempH, in_dimension = 1)
              #newExpo <- as.matrix(normalize_df_per_dim(tempH, 1))
              newSigs <- tempW * normFactor
              return(list(W = newSigs,
                          H = newExpo))
            })
            names(all_list) <- names(nmf_exp@WMatrix)

            # Set new matrices
            nmf_exp@WMatrix <- lapply(all_list, "[[" , "W")
            nmf_exp@HMatrix <- lapply(all_list, "[[" , "H")
            return(nmf_exp)
          }
)

#' @rdname regularizeH-methods
#' @aliases regularizeH,ANY,ANY-method
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 10)
#' regularizeH(nmf_exp)
#' }
setMethod("regularizeH",
          "ButchR_NMF",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            nmf_exp@WMatrix
            all_list <- lapply(seq_along(nmf_exp@WMatrix), function(k_ind){
              tempW <- nmf_exp@WMatrix[[k_ind]]
              tempH <- nmf_exp@HMatrix[[k_ind]]
              normFactor <- matrixStats::rowQuantiles(tempH, probs = 1)
              newExpo    <- tempH / normFactor
              newSigs    <- tempW * normFactor
              return(list(W = newSigs,
                          H = newExpo))
            })
            names(all_list) <- names(nmf_exp@WMatrix)

            # Set new matrices
            nmf_exp@WMatrix <- lapply(all_list, "[[" , "W")
            nmf_exp@HMatrix <- lapply(all_list, "[[" , "H")
            return(nmf_exp)
          }
)

#------------------------------------------------------------------------------#
#                       Signature specfific features                           #
#------------------------------------------------------------------------------#
# Signature specfific features

#' @rdname SignatureSpecificFeatures-methods
#' @aliases SignatureSpecificFeatures,ANY-method
#' @export
#'
#' @examples
#' \dontrun{
#' # For ButchR_NMF objects:
#' SignatureSpecificFeatures(nmf_exp)
#' lapply(SignatureSpecificFeatures(nmf_exp), function(x) sapply(x, length))
#' SignatureSpecificFeatures(nmf_exp, k = 3)
#' SignatureSpecificFeatures(nmf_exp, k = 3, return_all_features = TRUE)
#' }
setMethod("SignatureSpecificFeatures",
          "ButchR_NMF",
          function(x, k = NULL, return_all_features = FALSE, ...){
            # if no feature extraction was performed
            if (nrow(x@SignFeatures) == 0) {
              stop("No feature extraction has been performed, please run: \n",
                   "`compute_SignatureFeatures`")
            }

            # String of 0 and 1 to matrix
            bin_str_tu_mat <- function(binstr, return_all_features){
              names(binstr) <- rownames(x@SignFeatures)
              sig_feats <- as.matrix(do.call(rbind,
                                             lapply(strsplit(binstr, split = ""),
                                                    as.numeric)))
              if (return_all_features) {
                return(sig_feats)
              } else {
                sig_feats <- sig_feats[rowSums(sig_feats) == 1,,drop=FALSE]
                apply(sig_feats, 2, function(is_sig){
                  is_sig <- is_sig[is_sig == 1]
                  return(names(is_sig))
                })
              }
            }


            if(is.null(k)) {
              ssf <- lapply(x@SignFeatures, function(binstr) {
                bin_str_tu_mat(binstr, return_all_features)
              })
            } else {
              if (k == 2 ) {
                stop("Signature Specific Features extraction is not supported for K = 2")
              }

              idx <- as.character(x@OptKStats$rank_id[x@OptKStats$k == k])
              if (length(idx) == 0 ) {
                stop("No W matrix present for k = ", k,
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","))
              }
              ssf <- bin_str_tu_mat(x@SignFeatures[,idx], return_all_features)
            }
            return(ssf)
          }
)


#' @rdname compute_SignatureFeatures-methods
#' @aliases compute_SignatureFeatures,ANY-method
#'
#' #' @param x an object of class ButchR_NMF.
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 3,
#'                              method = "NMF",
#'                              n_initializations = 2,
#'                              extract_features = FALSE)
#' nmf_exp <- compute_SignatureFeatures(nmf_exp)
#' SignatureSpecificFeatures(nmf_exp, k = 3)
#' SignatureSpecificFeatures(nmf_exp, k = 3, return_all_features = TRUE)
#' }
setMethod("compute_SignatureFeatures",
          "ButchR_NMF",
          function(x){
            k <- x@OptKStats$k
            if (length(k) == 1 & 2 %in% k) {
              stop("Signature Specific Features extraction is not supported for K = 2")
            }

            #------------------------------------------------------------------#
            #             Compute signatures specific features                 #
            #------------------------------------------------------------------#
            SignFeatures <- lapply(x@WMatrix, function(W){
              if (ncol(W) == 2) {
                return(NULL)
              } else {
                rownames(W) <- x@input_matrix$rownames
                return(WcomputeFeatureStats(W))
              }
            })
            SignFeatures <- data.frame(do.call(cbind, SignFeatures),
                                       stringsAsFactors = FALSE)
            #------------------------------------------------------------------#
            #                     Return ButchR_NMF object                     #
            #------------------------------------------------------------------#
            x@SignFeatures <- SignFeatures
            return(x)
          }
)



