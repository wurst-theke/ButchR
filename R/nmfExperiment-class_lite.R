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
               FeatureStats              = "data.frame",
               #Factorization_ranks       = "data.frame",
               SignatureSpecificFeatures = "list" )

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


#==============================================================================#
#                                 Getter & Setter                              #
#==============================================================================#
#### H-Matrix (H-Matrix with smallest frobError)
# Getter
#setGeneric("HMatrix", function(x, k = NULL, ...) standardGeneric("HMatrix"))

#' Return a list of H-Matrices or an H-Matrix for the indicaded rank
#'
#' @param x an nmfExperiment_lite object
#' @param k numeric  - factorization rank
#'
#' @return list of H-Matrices or an H-Matrix for the indicaded rank
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
              idx <- x@OptKStats$rank_id[x@OptKStats$k == k]
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


#### W-Matrix (W-Matrix with smallest frobError)
# Getter
#setGeneric("WMatrix", function(x, k = NULL, ...) standardGeneric("WMatrix"))

#' Return a list of W-Matrices or a W-Matrix for the indicaded rank
#'
#' @param x an nmfExperiment_lite object
#' @param k numeric  - factorization rank
#'
#' @return list of W-Matrices or a W-Matrix for the indicaded rank
#' @export
#'
#' @examples
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
              idx <- x@OptKStats$rank_id[x@OptKStats$k == k]
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



#### FrobError
# Getter
#setGeneric("FrobError", function(x, ...) standardGeneric("FrobError"))

#' Return Frobenius Error from all initializations
#'
#' @param x an nmfExperiment_lite object
#'
#' @return data.frame with Frobenius error for all factorization ranks
#' @export
#'
#' @examples
#' FrobError(nmf_exp)
setMethod("FrobError", "nmfExperiment_lite", function(x, ...) x@FrobError)


#==============================================================================#
#                              Optimal K Statistics                            #
#==============================================================================#
#### Optimal K Statistics
# Getter
#setGeneric("OptKStats", function(x, ...) standardGeneric("OptKStats"))

#' Return optimal K Statistics
#'
#' @param x an nmfExperiment_lite object
#'
#' @return optimal K Statistics
#' @export
#'
#' @examples
#' OptKStats(nmf_exp)
setMethod("OptKStats", "nmfExperiment", function(x, ...) x@OptKStats)



#### Optimal K
# Getter
#setGeneric("OptK", function(x, ...) standardGeneric("OptK"))

#' Return optimal K
#'
#' @param x an nmfExperiment_lite object
#'
#' @return numeric - optimal K
#' @export
#'
#' @examples
#' OptK(nmf_exp)
setMethod("OptK", "nmfExperiment", function(x, ...) x@OptK)


















#' #### Feature Statistics
#' # Getter
#' setGeneric("FeatureStats", function(x, ...) standardGeneric("FeatureStats"))
#'
#' #' Feature Statistics getter
#' #'
#' #' @param nmfExperiment
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' setMethod("FeatureStats", "nmfExperiment", function(x, ...) x@FeatureStats)
#'
#' # Setter
#' setGeneric("setFeatureStats", function(nmfExperiment, FeatureStats)
#'   standardGeneric("setFeatureStats"))
#'
#' #' Feature Statistics setter
#' #'
#' #' @param nmfExperiment
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' setMethod("setFeatureStats", "nmfExperiment",
#'           function(nmfExperiment, FeatureStats) {
#'             if(nrow(nmfExperiment@FeatureStats) == 0) {
#'               nmfExperiment@FeatureStats <- FeatureStats
#'             } else {
#'               nmfExperiment@FeatureStats <-
#'                 cbind(nmfExperiment@FeatureStats, FeatureStats)
#'             }
#'             return(nmfExperiment)
#'           })
#'
#' #### Signature specfific features
#' # Getter
#' setGeneric("SignatureSpecificFeatures",
#'            function(x, ...) standardGeneric("SignatureSpecificFeatures"))
#'
#' #' Signature specfific features getter
#' #'
#' #' @param nmfExperiment
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' setMethod("SignatureSpecificFeatures",
#'           "nmfExperiment", function(x, ...) x@SignatureSpecificFeatures)
#'
#' # Setter
#' setGeneric("setSignatureSpecificFeatures",
#'            function(nmfExperiment, SignatureSpecificFeatures){
#'              standardGeneric("setSignatureSpecificFeatures")
#'            })
#'
#' #' Feature Statistics setter
#' #'
#' #' @param nmfExperiment
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' setMethod("setSignatureSpecificFeatures", "nmfExperiment",
#'           function(nmfExperiment, SignatureSpecificFeatures) {
#'             if(nrow(nmfExperiment@SignatureSpecificFeatures) == 0) {
#'               nmfExperiment@SignatureSpecificFeatures <- SignatureSpecificFeatures
#'             } else {
#'               nmfExperiment@SignatureSpecificFeatures <-
#'                 c(nmfExperiment@SignatureSpecificFeatures,
#'                   SignatureSpecificFeatures)
#'             }
#'             return(nmfExperiment)
#'           })
