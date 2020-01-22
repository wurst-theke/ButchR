# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#' NMF Experiment Class
#'
#' @slot HMatrixList list.
#' @slot HMatrix list.
#' @slot WMatrixList list.
#' @slot WMatrix list.
#' @slot FrobError DataFrame.
#' @slot OptKStats DataFrame.
#'
#' @return
#'
#' @import SummarizedExperiment
#' @export
#'
#' @examples
nmfExperiment <- setClass(
  Class = "nmfExperiment",
  contains = "SummarizedExperiment",
  representation = representation(HMatrixList = "list",
                                  WMatrixList = "list",
                                  FrobError = "DataFrame",
                                  OptKStats = "DataFrame",
                                  OptK = "numeric",
                                  FeatureStats = "DataFrame",
                                  SignatureSpecificFeatures = "list"))

#============================================================================#
#                                 Getter & Setter                            #
#============================================================================#
#### FrobError
# Getter
setGeneric("FrobError", function(x, ...) standardGeneric("FrobError"))

#' Frobenius Error getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("FrobError", "nmfExperiment", function(x, ...) x@FrobError)

# Setter
setGeneric("setFrobError", function(nmfExperiment, FrobError)
  standardGeneric("setFrobError"))

#' Frobenius Error setter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("setFrobError", "nmfExperiment", function(nmfExperiment, FrobError) {
  nmfExperiment@FrobError <- FrobError
  return(nmfExperiment)
})

#### H-Matrix List
# Getter
setGeneric("HMatrixList", function(x, k = NULL, ...)
  standardGeneric("HMatrixList"))

#' H-Matrix List getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("HMatrixList", "nmfExperiment", function(x, k = NULL, ...){
  if(is.null(k)) {
    x@HMatrixList
  } else {
    x@HMatrixList[[as.character(k)]]
  }
})

# Setter
setGeneric("setHMatrixList", function(nmfExperiment, HMatrixList)
  standardGeneric("setHMatrixList"))

#' H-Matrix List setter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("setHMatrixList", "nmfExperiment",
          function(nmfExperiment, HMatrixList){
            nmfExperiment@HMatrixList <- HMatrixList
            return(nmfExperiment)
          })

#### W-Matrix
# Getter
setGeneric("WMatrixList", function(x, k = NULL, ...)
  standardGeneric("WMatrixList"))

#' W-Matrix list getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("WMatrixList", "nmfExperiment", function(x, k = NULL, ...) {
  if(is.null(k)) {
    x@WMatrixList
  } else {
    x@WMatrixList[[as.character(k)]]
  }
})

# Setter
setGeneric("setWMatrixList", function(nmfExperiment, WMatrixList)
  standardGeneric("setWMatrixList"))

#' W-Matrix setter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("setWMatrixList", "nmfExperiment",
          function(nmfExperiment, WMatrixList) {
            nmfExperiment@WMatrixList <- WMatrixList
            return(nmfExperiment)
          })

#### H-Matrix (H-Matrix with smallest frobError)
# Getter
#' H-Matrix (H-Matrix with smallest frobError)
#'
#' Return a list of H-Matrices or an H-Matrix for the indicaded rank
#'
#' @param x an nmfExperiment or a nmfExperiment_lite object
#' @param k numeric  - factorization rank
#'
#' @return list of H-Matrices or an H-Matrix for the indicaded rank
#' @export
#' @docType methods
#' @rdname HMatrix-methods
#'
#' @examples
#' HMatrix(nmf_exp)
#' HMatrix(nmf_exp, k = 2)
setGeneric("HMatrix", function(x, k = NULL, ...)
  standardGeneric("HMatrix"))

#' H-Matrix getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("HMatrix", "nmfExperiment", function(x, k = NULL, ...) {
  i.min <- apply(x@FrobError, 2, which.min)
  if(is.null(k)) {
    H <- lapply(names(x@HMatrixList), function(k) {
      x@HMatrixList[[k]][[i.min[k]]]
    })
    names(H) <- names(x@HMatrixList)
  } else {
    k <- as.character(k)
    H <- x@HMatrixList[[k]][[i.min[k]]]
  }
  return(H)
})

#### W-Matrix (W-Matrix with smallest frobError)
# Getter
#' W-Matrix (W-Matrix with smallest frobError)
#'
#' Return a list of W-Matrices or a W-Matrix for the indicaded rank
#'
#' @param x an nmfExperiment or a nmfExperiment_lite object
#' @param k numeric  - factorization rank
#'
#' @return list of W-Matrices or a W-Matrix for the indicaded rank
#' @export
#' @docType methods
#' @rdname WMatrix-methods
#'
#' @examples
#' WMatrix(nmf_exp)
#' WMatrix(nmf_exp, k = 2)
setGeneric("WMatrix", function(x, k = NULL, ...) standardGeneric("WMatrix"))

#' W-Matrix getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("WMatrix", "nmfExperiment", function(x, k = NULL, ...) {
  i.min <- apply(x@FrobError, 2, which.min)
  if (is.null(k)) {
    W <- lapply(names(x@WMatrixList), function(k) {
      x@WMatrixList[[k]][[i.min[k]]]
    })
    names(W) <- names(x@WMatrixList)
  } else {
    k <- as.character(k)
    W <- x@WMatrixList[[k]][[i.min[k]]]
  }
  return(W)
})

#### Optimal K Statistics
# Getter
setGeneric("OptKStats", function(x, ...) standardGeneric("OptKStats"))

#' Optimal K Statistics getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("OptKStats", "nmfExperiment", function(x, ...) x@OptKStats)

# Setter
setGeneric("setOptKStats", function(nmfExperiment, OptKStats)
  standardGeneric("setOptKStats"))

#' Optimal K Statistics setter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("setOptKStats", "nmfExperiment",
          function(nmfExperiment, OptKStats) {
            nmfExperiment@OptKStats <- OptKStats
            return(nmfExperiment)
          })

#### Optimal K
# Getter
setGeneric("OptK", function(x, ...) standardGeneric("OptK"))

#' Optimal K
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("OptK", "nmfExperiment", function(x, ...) x@OptK)

# Setter
setGeneric("setOptK", function(nmfExperiment, OptK) standardGeneric("setOptK"))

#' Optimal K setter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("setOptK", "nmfExperiment", function(nmfExperiment, OptK) {
  nmfExperiment@OptK <- OptK
  return(nmfExperiment)
})

#### Feature Statistics
# Getter
setGeneric("FeatureStats", function(x, ...) standardGeneric("FeatureStats"))

#' Feature Statistics getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("FeatureStats", "nmfExperiment", function(x, ...) x@FeatureStats)

# Setter
setGeneric("setFeatureStats", function(nmfExperiment, FeatureStats)
  standardGeneric("setFeatureStats"))

#' Feature Statistics setter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("setFeatureStats", "nmfExperiment",
          function(nmfExperiment, FeatureStats) {
            if(nrow(nmfExperiment@FeatureStats) == 0) {
              nmfExperiment@FeatureStats <- FeatureStats
            } else {
              nmfExperiment@FeatureStats <-
                cbind(nmfExperiment@FeatureStats, FeatureStats)
            }
            return(nmfExperiment)
          })

#### Signature specfific features
# Getter
setGeneric("SignatureSpecificFeatures",
           function(x, ...) standardGeneric("SignatureSpecificFeatures"))

#' Signature specfific features getter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("SignatureSpecificFeatures",
          "nmfExperiment", function(x, ...) x@SignatureSpecificFeatures)

# Setter
setGeneric("setSignatureSpecificFeatures",
           function(nmfExperiment, SignatureSpecificFeatures){
             standardGeneric("setSignatureSpecificFeatures")
           })

#' Feature Statistics setter
#'
#' @param nmfExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("setSignatureSpecificFeatures", "nmfExperiment",
          function(nmfExperiment, SignatureSpecificFeatures) {
  if(nrow(nmfExperiment@SignatureSpecificFeatures) == 0) {
    nmfExperiment@SignatureSpecificFeatures <- SignatureSpecificFeatures
  } else {
    nmfExperiment@SignatureSpecificFeatures <-
      c(nmfExperiment@SignatureSpecificFeatures,
        SignatureSpecificFeatures)
  }
  return(nmfExperiment)
})



#==============================================================================#
#                               NMF Normalization                              #
#==============================================================================#
#' Normalize the signatures matrix (W)
#'
#' After column normalization of the matrix W, the inverse factors are
#' mutiplied with the rows of H in order to keep the matrix product W*H
#' constant.
#'
#' @param nmf_exp an nmfExperiment or a nmfExperiment_lite object
#'
#' @return an nmfExperiment or a nmfExperiment_lite object normalized by W
#' @export
#' @docType methods
#' @rdname normalizeW-methods
#'
#' @examples
#' normalizeW(nmf_exp)
setGeneric("normalizeW", function(nmf_exp, ...) standardGeneric("normalizeW"))

#' @rdname normalizeW-methods
#' @aliases normalizeW,ANY,ANY-method
#'
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
setMethod("normalizeW", "nmfExperiment",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            all_list <- lapply(seq_along(WMatrixList(nmf_exp)), function(k_ind){
              k_list <-
                lapply(seq_along(WMatrixList(nmf_exp)[[k_ind]]), function(init_ind){
                  tempW <- WMatrixList(nmf_exp)[[k_ind]][[init_ind]]
                  tempH <- HMatrixList(nmf_exp)[[k_ind]][[init_ind]]
                  normFactor <- colSums(tempW)
                  # catch errors associated with NaNs in W or H
                  if (any(is.nan(normFactor))){
                    return(list(W = tempW,
                                H = tempH))
                  }else{
                    newSigs <- as.matrix(normalize_df_per_dim(tempW, 2))
                    newExpo <- tempH * normFactor
                    #newV <- newSigs %*% newExpo
                    #oldV <- tempW %*% tempH
                    return(list(W = newSigs,
                                H = newExpo))
                  }
                })
              names(k_list) <- names(WMatrixList(nmf_exp)[[k_ind]])
              return(k_list)
            })
            names(all_list) <- names(WMatrixList(nmf_exp))
            thisWMatrixList <- lapply(all_list, function(current_k_list){
              kWMatrixList <- lapply(current_k_list, function(current_entry){
                return(current_entry$W)
              })
            })
            nmf_exp <- setWMatrixList(nmf_exp, thisWMatrixList)
            thisHMatrixList <- lapply(all_list, function(current_k_list){
              kHMatrixList <- lapply(current_k_list, function(current_entry){
                return(current_entry$H)
              })
            })
            nmf_exp <- setHMatrixList(nmf_exp, thisHMatrixList)
            return(nmf_exp)
          }
)




#' Normalize the signatures matrix (H)
#'
#' After row normalization of the matrix H, the inverse factors are
#' mutiplied with the columns of W in order to keep the matrix product W*H
#' constant.
#'
#' @param nmf_exp an nmfExperiment or a nmfExperiment_lite object
#'
#' @return an nmfExperiment or a nmfExperiment_lite object normalized by W
#' @export
#' @docType methods
#' @rdname normalizeH-methods
#'
#' @examples
#' normalizeH(nmf_exp)
setGeneric("normalizeH", function(nmf_exp, ...) standardGeneric("normalizeH"))

#' @rdname normalizeH-methods
#' @aliases normalizeH,ANY,ANY-method
#'
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
setMethod("normalizeH", "nmfExperiment",
          function(nmf_exp){
            # account for WMatrixList and HMatrixList
            all_list <- lapply(seq_along(WMatrixList(nmf_exp)), function(k_ind){
              k_list <-
                lapply(seq_along(WMatrixList(nmf_exp)[[k_ind]]), function(init_ind){
                  tempW <- WMatrixList(nmf_exp)[[k_ind]][[init_ind]]
                  tempH <- HMatrixList(nmf_exp)[[k_ind]][[init_ind]]
                  normFactor <- rowSums(tempH)
                  newExpo <- as.matrix(normalize_df_per_dim(tempH, 1))
                  newSigs <- tempW * normFactor
                  return(list(W = newSigs,
                              H = newExpo))
                })
              names(k_list) <- names(WMatrixList(nmf_exp)[[k_ind]])
              return(k_list)
            })
            names(all_list) <- names(WMatrixList(nmf_exp))
            thisWMatrixList <- lapply(all_list, function(current_k_list){
              kWMatrixList <- lapply(current_k_list, function(current_entry){
                return(current_entry$W)
              })
            })
            nmf_exp <- setWMatrixList(nmf_exp, thisWMatrixList)
            thisHMatrixList <- lapply(all_list, function(current_k_list){
              kHMatrixList <- lapply(current_k_list, function(current_entry){
                return(current_entry$H)
              })
            })
            nmf_exp <- setHMatrixList(nmf_exp, thisHMatrixList)
            return(nmf_exp)
          }
)


#' Regularize the signatures matrix (H)
#'
#' After row regularization of the matrix H, the inverse factors are
#' mutiplied with the columns of W in order to keep the matrix product W*H
#' constant.
#'
#' @param nmf_exp
#'
#' @return A data structure of type nmfExperiment
#'
#' @export
#'
#' @examples
#'  NULL
#'
regularizeH <- function(nmf_exp){
  # account for WMatrixList and HMatrixList
  all_list <- lapply(seq_along(WMatrixList(nmf_exp)), function(k_ind){
    k_list <-
      lapply(seq_along(WMatrixList(nmf_exp)[[k_ind]]), function(init_ind){
        tempW <- WMatrixList(nmf_exp)[[k_ind]][[init_ind]]
        tempH <- HMatrixList(nmf_exp)[[k_ind]][[init_ind]]
        normFactor <- rowMax(tempH)
        newExpo <- tempH / normFactor
        newSigs <- tempW * normFactor
        return(list(W = newSigs,
                    H = newExpo))
      })
    names(k_list) <- names(WMatrixList(nmf_exp)[[k_ind]])
    return(k_list)
  })
  names(all_list) <- names(WMatrixList(nmf_exp))
  thisWMatrixList <- lapply(all_list, function(current_k_list){
    kWMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$W)
    })
  })
  nmf_exp <- setWMatrixList(nmf_exp, thisWMatrixList)
  thisHMatrixList <- lapply(all_list, function(current_k_list){
    kHMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$H)
    })
  })
  nmf_exp <- setHMatrixList(nmf_exp, thisHMatrixList)
  return(nmf_exp)
}
