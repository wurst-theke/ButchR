#' @include generics.R
NULL

#------------------------------------------------------------------------------#
#                  Join NMF tensorflow Wrapper - CLASS - FUNCTION              #
#------------------------------------------------------------------------------#

#' Join NMF Class
#'
#' @slot HMatrix list.
#' @slot WMatrix_vs list.
#' @slot FrobError DataFrame.
#' @slot OptKStats DataFrame.
#'
#' @return
#'
#' @export
#'
#' @examples
join_NMF <- setClass(
  Class = "join_NMF",
  slots = list(input_data   = "list",
               HMatrix      = "list",
               WMatrix_vs   = "list",
               FrobError    = "data.frame",
               OptKStats    = "data.frame",
               OptK         = "numeric",
               FeatureStats              = "data.frame",
               SignatureSpecificFeatures = "list")
)

setMethod("show",
          "join_NMF",
          function(object) {
            #cat("class: join NMF object \n")
            #cat("Best factorization index: ", object@best_factorization_idx, "\n")
            cat("class: join NMF object \n")
            x <- apply(object@input_data$dim, 1, paste, collapse = ",")
            x <- gsub(" ", "", x)
            cat("Original matrices dimensions: ", x, "\n")
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
#                                 Getter & Setter                              #
#------------------------------------------------------------------------------#
#### H-Matrix (H-Matrix with smallest frobError)
#' Return a list of H-Matrices or an H-Matrix for the indicaded rank
#'
#' @param x a join_NMF object
#' @param k numeric - factorization rank
#'
#' @return list of H-Matrices or an H-Matrix for the indicaded rank
#' @export
#'
#' @examples
#' HMatrix(jnmf_exp)
#' HMatrix(jnmf_exp, k = 2)
setMethod("HMatrix",
          "join_NMF",
          function(x, k = NULL, ...) {
            if(is.null(k)) {
              H <- x@HMatrix
              if (!is.null(x@input_data$colnames)) {
                H <- lapply(H, function(hi) {
                  colnames(hi) <- x@input_data$colnames
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
              if (!is.null(x@input_data$colnames)) {
                colnames(H) <- x@input_data$colnames
              }
            }
            return(H)
          }
)

#### W-Matrix (W-Matrix with smallest frobError)
# Getter
#' Return a list of W-Matrices or a W-Matrix for the indicaded rank
#'
#' @param x an join_NMF object
#' @param k numeric  - factorization rank
#'
#' @return list of W-Matrices or a W-Matrix for the indicaded rank
#' @export
#'
#' @examples
#' WMatrix(jnmf_exp)
#' WMatrix(jnmf_exp, k = 2)
setMethod("WMatrix",
          "join_NMF",
          function(x, k = NULL, ...) {
            if(is.null(k)) {
              W <- lapply(x@WMatrix_vs, function(wi) {
                view_ids <- setNames(names(wi), names(wi))
                lapply(view_ids, function(view_id){
                  rownames(wi[[view_id]]) <- x@input_data$rownames[[view_id]]
                  wi[[view_id]]
                })

              })
            } else {
              idx <- x@OptKStats$rank_id[x@OptKStats$k == k]
              if (length(idx) == 0 ) {
                stop("No W matrix present for k = ", k,
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","))
              }
              W <- x@WMatrix_vs[[idx]]
              view_ids <- setNames(names(W), names(W))
              W <- lapply(view_ids, function(view_id){
                rownames(W[[view_id]]) <- x@input_data$rownames[[view_id]]
                W[[view_id]]
              })
            }
            return(W)
          }
)



#' Return Frobenius Error from all initializations
#'
#' @param x an nmfExperiment_lite object
#'
#' @return data.frame with Frobenius error for all factorization ranks
#' @export
#'
#' @examples
#' FrobError(jnmf_exp)
setMethod("FrobError", "join_NMF", function(x, ...) x@FrobError)


