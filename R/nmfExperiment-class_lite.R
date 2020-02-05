#' @include generics.R
NULL

# Copyright © 2015-2020  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#' NMF Experiment Class lite
#'
#' @slot HMatrix list.
#' @slot WMatrix list.
#' @slot FrobError DataFrame.
#' @slot OptKStats DataFrame.
#'
#' @return An object of NMF Experiment Class lite
#'
#' @export
#'
#' @examples
nmfExperiment_lite <- setClass(
  Class = "nmfExperiment_lite",
  slots = list(input_matrix = "list",
               HMatrix      = "list",
               WMatrix      = "list",
               FrobError    = "data.frame",
               OptKStats    = "data.frame",
               OptK         = "numeric",
               SignFeatures = "data.frame" )

)

setMethod("show",
          "nmfExperiment_lite",
          function(object) {
            cat("class: NMF object lite \n")
            #cat("Best factorization index: ", object@best_factorization_idx, "\n")
            cat("Original matrix dimension: ", object@input_matrix$dim, "\n")
            cat("Factorization performed for ranks: ", object@OptKStats$k, "\n")
            if (length(object@OptK) == 0) {
              OptK <- "Please select manualy\n"
            } else {
              OptK <- object@OptK
            }

            cat("Optimal K based on factorization metrics: ", OptK, "\n")
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
          "nmfExperiment_lite",
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
#' # For nmfExperiment_lite objects:
#' WMatrix(nmf_exp)
#' WMatrix(nmf_exp, k = 2)
setMethod("WMatrix",
          "nmfExperiment_lite",
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
#' @aliases SignatureSpecificFeatures,ANY-method
#' @export
#'
setMethod("FrobError", "nmfExperiment_lite", function(x, ...) x@FrobError)


#### Optimal K Statistics
# Getter
#' Return optimal factorization rank (K) Statistics
#'
#' @param x an nmfExperiment_lite object
#'
#' @return optimal K Statistics
#' @export
#'
#' @examples
#' OptKStats(nmf_exp)
setMethod("OptKStats", "nmfExperiment_lite", function(x, ...) x@OptKStats)

#### Optimal K
# Getter
#' Return optimal K
#'
#' @param x an nmfExperiment_lite object
#'
#' @return numeric - optimal K
#' @export
#'
#' @examples
#' OptK(nmf_exp)
setMethod("OptK", "nmfExperiment_lite", function(x, ...) x@OptK)

#------------------------------------------------------------------------------#
#                               NMF Normalization                              #
#------------------------------------------------------------------------------#
#' @rdname normalizeW-methods
#' @aliases normalizeW,ANY,ANY-method
#'
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
setMethod("normalizeW",
          "nmfExperiment_lite",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            #nmf_exp@WMatrix
            all_list <- lapply(seq_along(nmf_exp@WMatrix), function(k_ind){
              tempW <- nmf_exp@WMatrix[[k_ind]]
              tempH <- nmf_exp@HMatrix[[k_ind]]
              normFactor <- colSums(tempW)
              newSigs <- as.matrix(YAPSA::normalize_df_per_dim(tempW, 2))
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
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
setMethod("normalizeH",
          "nmfExperiment_lite",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            #nmf_exp@WMatrix
            all_list <- lapply(seq_along(nmf_exp@WMatrix), function(k_ind){
              tempW <- nmf_exp@WMatrix[[k_ind]]
              tempH <- nmf_exp@HMatrix[[k_ind]]
              normFactor <- rowSums(tempH)
              newExpo <- as.matrix(YAPSA::normalize_df_per_dim(tempH, 1))
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
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
setMethod("regularizeH",
          "nmfExperiment_lite",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            nmf_exp@WMatrix
            all_list <- lapply(seq_along(nmf_exp@WMatrix), function(k_ind){
              tempW <- nmf_exp@WMatrix[[k_ind]]
              tempH <- nmf_exp@HMatrix[[k_ind]]
              normFactor <- rowMax(tempH)
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
#' # For nmfExperiment_lite objects:
#' SignatureSpecificFeatures(nmf_exp)
#' lapply(SignatureSpecificFeatures(nmf_exp), function(x) sapply(x, length))
#' SignatureSpecificFeatures(nmf_exp, k = 3)
#' SignatureSpecificFeatures(nmf_exp, k = 3, return_all_features = TRUE)
setMethod("SignatureSpecificFeatures",
          "nmfExperiment_lite",
          function(x, k = NULL, return_all_features = FALSE, ...){

            # String of 0 and 1 to matrix
            bin_str_tu_mat <- function(binstr, return_all_features){
              names(binstr) <- rownames(x@SignFeatures)
              sig_feats <- do.call(rbind, lapply(strsplit(binstr, split = ""), as.numeric))
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
