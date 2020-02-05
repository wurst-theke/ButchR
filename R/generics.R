#### FrobError
#' Frobenius error for all factorization ranks
#'
#' Returns a data.frame with factorization ranks in the columns,
#' and the Frobenius error for every NMF initialization in the rows.
#'
#' @param x an R object of type nmfExperiment, nmfExperiment_lite, join_NMF or
#' integrative_NMF.
#'
#' @return data.frame with Frobenius errors
#' @docType methods
#' @rdname FrobError-methods
#' @exportMethod FrobError
#' @examples
#' FrobError(nmf_exp)
setGeneric("FrobError", function(x, ...) standardGeneric("FrobError"))

# Setter
setGeneric("setFrobError", function(nmfExperiment, FrobError)
  standardGeneric("setFrobError"))

#### H-Matrix List
# Getter
setGeneric("HMatrixList", function(x, k = NULL, ...)
  standardGeneric("HMatrixList"))

# Setter
setGeneric("setHMatrixList", function(nmfExperiment, HMatrixList)
  standardGeneric("setHMatrixList"))


#### W-Matrix
# Getter
setGeneric("WMatrixList", function(x, k = NULL, ...)
  standardGeneric("WMatrixList"))

# Setter
setGeneric("setWMatrixList", function(nmfExperiment, WMatrixList)
  standardGeneric("setWMatrixList"))


#' H-Matrix (H-Matrix with smallest frobError)
#'
#' Return a list of H-Matrices or an H-Matrix for the indicaded rank
#'
#' @param x an R object of type nmfExperiment, nmfExperiment_lite, join_NMF or
#' integrative_NMF.
#' @param k numeric - factorization rank.
#'
#' @return list of H-Matrices or an H-Matrix for the indicaded rank.
#' @docType methods
#' @rdname HMatrix-methods
#' @exportMethod HMatrix
#' @examples
#' HMatrix(nmf_exp)
#' HMatrix(nmf_exp, k = 2)
setGeneric("HMatrix", function(x, k = NULL, ...)
  standardGeneric("HMatrix"))


#' W-Matrix (W-Matrix with smallest frobError)
#'
#' Return a list of W-Matrices or a W-Matrix for the indicaded rank
#'
#' @param x an R object of type nmfExperiment, nmfExperiment_lite, join_NMF or
#' integrative_NMF.
#' @param k numeric - factorization rank.
#'
#' @return list of W-Matrices or a W-Matrix for the indicaded rank.
#' @docType methods
#' @rdname WMatrix-methods
#' @exportMethod WMatrix
#' @examples
#' WMatrix(nmf_exp)
#' WMatrix(nmf_exp, k = 2)
setGeneric("WMatrix", function(x, k = NULL, ...) standardGeneric("WMatrix"))



#### Optimal K Statistics
# Getter
setGeneric("OptKStats", function(x, ...) standardGeneric("OptKStats"))

# Setter
setGeneric("setOptKStats", function(nmfExperiment, OptKStats)
  standardGeneric("setOptKStats"))



#### Optimal K
# Getter
setGeneric("OptK", function(x, ...) standardGeneric("OptK"))

# Setter
setGeneric("setOptK", function(nmfExperiment, OptK) standardGeneric("setOptK"))




#### Feature Statistics
# Getter
setGeneric("FeatureStats", function(x, ...) standardGeneric("FeatureStats"))

# Setter
setGeneric("setFeatureStats", function(nmfExperiment, FeatureStats)
  standardGeneric("setFeatureStats"))


#### Signature specfific features
#' Signature Specific Features
#'
#' Returns the list of signatures specific features
#' for all factorization ranks or for the indicaded rank,
#' if return_all_features = TRUE
#' returns a binary matrix for every factorization rank,
#' with features in the rows and samples in the columns,
#' in which 1 means that the features is contributing to the signature,
#' and 0 it does not.
#' The extraction of Signature Specific Features is not supported for k = 2
#'
#' @param x an nmfExperiment or a nmfExperiment_lite object
#' @param k numeric  - factorization rank
#'
#' @return list of signature specific fatures or binary matrices for all features
#' @docType methods
#' @rdname SignatureSpecificFeatures-methods
#' @exportMethod SignatureSpecificFeatures
#' @examples
#' SignatureSpecificFeatures(nmf_exp)
#' SignatureSpecificFeatures(nmf_exp, k = 2)
#' SignatureSpecificFeatures(nmf_exp, k = 2, return_all_features = TRUE)
setGeneric("SignatureSpecificFeatures",
           function(x, k = NULL, return_all_features = FALSE, ...)
             standardGeneric("SignatureSpecificFeatures"))

# Setter
setGeneric("setSignatureSpecificFeatures",
           function(nmfExperiment, SignatureSpecificFeatures){
             standardGeneric("setSignatureSpecificFeatures")
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

#' Regularize the signatures matrix (H)
#'
#' After row regularization of the matrix H, the inverse factors are
#' mutiplied with the columns of W in order to keep the matrix product W*H
#' constant.
#'
#' @param nmf_exp
#'
#' @return an nmfExperiment or a nmfExperiment_lite object normalized by W
#'
#' @export
#' @docType methods
#' @rdname regularizeH-methods
#'
#' @examples
#' regularizeH(nmf_exp)
setGeneric("regularizeH", function(nmf_exp, ...) standardGeneric("regularizeH"))