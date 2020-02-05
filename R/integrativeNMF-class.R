#' @include nmfExperiment-class_lite.R
NULL

#------------------------------------------------------------------------------#
#             Integrative NMF tensorflow Wrapper - CLASS - FUNCTION            #
#------------------------------------------------------------------------------#

#' Integrative NMF Class
#'
#' @slot HMatrix list.
#' @slot HMatrix_vs list.
#' @slot WMatrix_vs list.
#' @slot FrobError DataFrame.
#' @slot OptKStats DataFrame.
#'
#' @return
#'
#' @export
#'
#' @examples
integrative_NMF <- setClass(
  Class = "integrative_NMF",
  slots = list(input_data   = "list",
               HMatrix      = "list",
               HMatrix_vs   = "list",
               WMatrix_vs   = "list",
               FrobError    = "data.frame",
               OptKStats    = "data.frame",
               OptK         = "numeric",
               SignFeatures = "list")
)

setMethod("show",
          "integrative_NMF",
          function(object) {
            #cat("class: integrative NMF object \n")
            #cat("Best factorization index: ", object@best_factorization_idx, "\n")
            cat("class: integrative NMF object \n")
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
            rparams <- object@input_data$run_params
            cat("Running parameters: \n")
            cat("method = integrative NMF \n")
            cat("lamb = ", rparams$lamb, " \n")
            cat("n_initializations = ", rparams$n_initializations, " \n")
            cat("iterations = ", rparams$iterations, " \n")
            cat("stop_threshold = ", rparams$stop_threshold, " \n")
            cat("extract_features = ", rparams$extract_features, " \n")
          }
)


#------------------------------------------------------------------------------#
#                               H and W matrices                               #
#------------------------------------------------------------------------------#
#### H-Matrix (H-Matrix with smallest frobError)
#' @rdname HMatrix-methods
#' @aliases HMatrix,ANY-method
#'
#' @param view_id for integrative NMF; character vector with views
#'  from which to extract H matrices
#' @param type for integrative NMF; type of H matrix to extract, could be:
#' \itemize{
#' \item shared - shared H matrix
#' \item viewspec - view specific H matrix
#' \item total - sum of shared H matrix and view specific H matrix.
#' \item all - shared H matrix and view specific H matrices.
#' }
#' @export
#'
#' @examples
#' # For integrative_NMF objects:
#' # extract H matrices for all factorization ranks
#' HMatrix(inmf_exp, type = "shared")
#' HMatrix(inmf_exp, type = "viewspec")
#' HMatrix(inmf_exp, type = "total")
#' HMatrix(inmf_exp, type = "all")
#' # extract H matrices only for selected rank
#' HMatrix(inmf_exp, k = 2, type = "shared")
#' HMatrix(inmf_exp, k = 2, type = "viewspec")
#' HMatrix(inmf_exp, k = 2, type = "total")
#' HMatrix(inmf_exp, k = 2, type = "all")
#' # extract H matrices only for selected view and rank
#' HMatrix(inmf_exp, k = 2, view_id = "atac", type = "viewspec")
#' HMatrix(inmf_exp, k = 2, view_id = "atac", type = "total")
#'
setMethod("HMatrix",
          "integrative_NMF",
          function(x, k = NULL, view_id = NULL, type, ...) {
            # Check if view id is indeed one of the views
            if (is.null(view_id)) {
              view_id <- as.character(x@input_data$dim$view_ids)
              view_id <- setNames(view_id, view_id)
            } else if (all(view_id %in% x@input_data$dim$view_ids)) {
              view_id <- setNames(view_id, view_id)
            } else {
              view_id <- view_id[!view_id %in% x@input_data$dim$view_ids]
              stop("View: ", paste0(view_id, collapse = ","),
                   "\nis not present, please select from views = ",
                   paste0(x@input_data$dim$view_ids, collapse = ","))
            }
            # Only selected views and assign colnames
            Hshared <- lapply(x@HMatrix, function(hi) {
              colnames(hi) <- x@input_data$colnames
              hi
            })

            HMatrix_vs <- lapply(x@HMatrix_vs, function(hi) {
              lapply(view_id, function(view_idi){
                colnames(hi[[view_idi]]) <- x@input_data$colnames
                hi[[view_idi]]
              })
            })


            if (type == "shared") {
              Hfinal <- Hshared
            } else if (type == "viewspec") {
              Hfinal <- HMatrix_vs
            } else if (type == "total") {
              rank_ids <- setNames(names(HMatrix_vs), names(HMatrix_vs))
              Hfinal <- lapply(rank_ids, function(rank_id){
                lapply(HMatrix_vs[[rank_id]], function(hvs_rank){
                  Hshared[[rank_id]] + hvs_rank
                })
              })

            } else if (type == "all") {
              rank_ids <- setNames(names(HMatrix_vs), names(HMatrix_vs))
              Hfinal <- lapply(rank_ids, function(rank_id){
                c(list(shared = Hshared[[rank_id]]), HMatrix_vs[[rank_id]])
              })

            } else {
              stop("Please select the type of H matrices to retrieve from:",
                   "\n type = 'shared' 'viewspec' 'total' 'all'")
            }

            # Selected ranks
            if(!is.null(k)) {
              idx <- x@OptKStats$rank_id[x@OptKStats$k == k]
              if (length(idx) == 0 ) {
                stop("No H matrix present for k = ", k,
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","))
              }
              Hfinal <- Hfinal[[idx]]
            }
            # return only one matrix id list is equal to 1
            if (length(Hfinal) == 1) {
              Hfinal <- Hfinal[[1]]
            }
            return(Hfinal)
          }
)


# W-Matrix (W-Matrix with smallest frobError)
#' @rdname WMatrix-methods
#' @aliases WMatrix,ANY-method
#' @export
#'
#' @examples
#' # For integrative_NMF objects:
#' WMatrix(inmf_exp)
#' WMatrix(inmf_exp, k = 2)
#' lapply(WMatrix(inmf_exp, k = 2), head)
#' WMatrix(inmf_exp, k = 2, view_id = "atac")
setMethod("WMatrix",
          "integrative_NMF",
          function(x, k = NULL, view_id = NULL, ...) {
            # Check if view id is indeed one of the views
            if (is.null(view_id)) {
              view_id <- as.character(x@input_data$dim$view_ids)
              view_id <- setNames(view_id, view_id)
            } else if (all(view_id %in% x@input_data$dim$view_ids)) {
              view_id <- setNames(view_id, view_id)
            } else {
              view_id <- view_id[!view_id %in% x@input_data$dim$view_ids]
              stop("View: ", paste0(view_id, collapse = ","),
                   "\nis not present, please select from views = ",
                   paste0(x@input_data$dim$view_ids, collapse = ","))
            }

            # Only selected views
            WMatrix_vs <- lapply(x@WMatrix_vs, function(wi) {
              lapply(view_id, function(view_idi){
                rownames(wi[[view_idi]]) <- x@input_data$rownames[[view_idi]]
                wi[[view_idi]]
              })
            })

            if(is.null(k)) {
              W <- WMatrix_vs
            } else {
              idx <- x@OptKStats$rank_id[x@OptKStats$k == k]
              if (length(idx) == 0 ) {
                stop("No W matrix present for k = ", k,
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","))
              }
              W <- WMatrix_vs[[idx]]
            }
            # return only one matrix id list is equal to 1
            if (length(W) == 1) {
              W <- W[[1]]
            }
            return(W)
          }
)

# Return Frobenius Error from all initializations

#' @rdname FrobError-methods
#' @aliases SignatureSpecificFeatures,ANY-method
#' @export
#'
setMethod("FrobError", "integrative_NMF", function(x, ...) x@FrobError)


#------------------------------------------------------------------------------#
#                       Signature specfific features                           #
#------------------------------------------------------------------------------#
# Returns Signature specfific features
# For join NMF view_id is the name of one of original matrix to retrieve features from

#' @rdname SignatureSpecificFeatures-methods
#' @aliases SignatureSpecificFeatures,ANY-method
#'
#' @param view_id character vector with views from which sigature specific
#' features will be extracted.
#' @export
#'
#' @examples
#' # For integrative_NMF objects:
#' SignatureSpecificFeatures(inmf_exp)
#' lapply(SignatureSpecificFeatures(inmf_exp), function(view) sapply(view, function(x) sapply(x, length)))
#' lapply(SignatureSpecificFeatures(inmf_exp, k = 3), function(view) sapply(view, length))
#' SignatureSpecificFeatures(inmf_exp, k = 3, return_all_features = TRUE)
#' SignatureSpecificFeatures(inmf_exp, k = 3, return_all_features = TRUE, view_id = "atac")
#' SignatureSpecificFeatures(inmf_exp, return_all_features = TRUE, view_id = "atac")
setMethod("SignatureSpecificFeatures",
          "integrative_NMF",
          function(x, k = NULL, return_all_features = FALSE, view_id = NULL, ...){
            # Check if view id is indeed one of the views
            if (is.null(view_id)) {
              view_id <- setNames(names(x@SignFeatures), names(x@SignFeatures))
            } else if (all(view_id %in% names(x@SignFeatures))) {
              view_id <- setNames(view_id, view_id)
            } else {
              view_id <- view_id[!view_id %in% names(x@SignFeatures)]
              stop("View: ", paste0(view_id, collapse = ","),
                   "\nis not present, please select from views = ",
                   paste0(names(x@SignFeatures), collapse = ","))
            }
            # Only selected views
            SignFeatures <- x@SignFeatures[view_id]


            # String of 0 and 1 to matrix
            bin_str_tu_mat <- function(binstr, return_all_features, feature_ids){
              names(binstr) <- feature_ids
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


            # Extract for all ranks
            if(is.null(k)) {
              ssf <- lapply(SignFeatures, function(viewSignF){
                #print(head(viewSignF))
                lapply(viewSignF, function(binstr) {
                  bin_str_tu_mat(binstr, return_all_features, rownames(viewSignF))
                })
              })
              # Extract for selected rank
            } else {
              if (k == 2 ) {
                stop("Signature Specific Features extraction is not supported for K = 2")
              }
              idx <- as.character(x@OptKStats$rank_id[x@OptKStats$k == k])
              if (length(idx) == 0 ) {
                stop("No W matrix present for k = ", k,
                     "\nPlease select from ranks = ", paste0(x@OptKStats$k, collapse = ","))
              }

              ssf <- lapply(SignFeatures, function(viewSignF){
                bin_str_tu_mat(viewSignF[,idx], return_all_features, rownames(viewSignF))
              })
            }
            return(ssf)
          }
)


