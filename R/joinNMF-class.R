#' @include nmfExperiment-class_lite.R
NULL

#------------------------------------------------------------------------------#
#                  Join NMF tensorflow Wrapper - CLASS - FUNCTION              #
#------------------------------------------------------------------------------#

#' Join NMF Class
#' @slot input_data lsit of matrices
#' @slot HMatrix list.
#' @slot WMatrix_vs list.
#' @slot FrobError DataFrame.
#' @slot OptKStats DataFrame.
#' @slot OptK numeric
#' @slot SignFeatures DataFrame or list
#'
#' @return An object of class ButchR_joinNMF.
#' @export
#'
#' @examples
#' \dontrun{
#' ButchR_joinNMF(input_data   = input_data,
#'                HMatrix      = shared_HMatrix_list,
#'                WMatrix_vs   = view_specific_WMatrix_list,
#'                FrobError    = frob_errors,
#'                OptKStats    = OptKStats,
#'                OptK         = OptK,
#'                SignFeatures = SignFeatures)
#' }
ButchR_joinNMF <- setClass(
  Class = "ButchR_joinNMF",
  slots = list(input_data   = "list",
               HMatrix      = "list",
               WMatrix_vs   = "list",
               FrobError    = "data.frame",
               OptKStats    = "data.frame",
               OptK         = "numeric",
               SignFeatures = "list")
)

setMethod("show",
          "ButchR_joinNMF",
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
            rparams <- object@input_data$run_params
            cat("Running parameters: \n")
            cat("method = join NMF \n")
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
#' @export
#'
#' @examples
#' \dontrun{
#' # For ButchR_joinNMF objects:
#' HMatrix(jnmf_exp)
#' HMatrix(jnmf_exp, k = 2)
#' lapply(HMatrix(jnmf_exp, k = 2), head)
#' HMatrix(jnmf_exp, k = 2, view_id = "atac")
#' }
setMethod("HMatrix",
          "ButchR_joinNMF",
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

# W-Matrix (W-Matrix with smallest frobError)
#' @rdname WMatrix-methods
#' @aliases WMatrix,ANY-method
#'
#' @param view_id character vector with views from which to extract W matrices
#' @export
#'
#' @examples
#' \dontrun{
#' # For ButchR_joinNMF objects:
#' WMatrix(jnmf_exp)
#' WMatrix(jnmf_exp, k = 2)
#' lapply(WMatrix(jnmf_exp, k = 2), head)
#' WMatrix(jnmf_exp, k = 2, view_id = "atac")
#' }
setMethod("WMatrix",
          "ButchR_joinNMF",
          function(x, k = NULL, view_id = NULL, ...) {
            # Check if view id is indeed one of the views
            if (is.null(view_id)) {
              view_id <- as.character(x@input_data$dim$view_ids)
              view_id <- stats::setNames(view_id, view_id)
            } else if (all(view_id %in% x@input_data$dim$view_ids)) {
              view_id <- stats::setNames(view_id, view_id)
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
#' @aliases FrobError,ANY-method
#' @export
#'
setMethod("FrobError", "ButchR_joinNMF", function(x, ...) x@FrobError)


#------------------------------------------------------------------------------#
#                       Signature specfific features                           #
#------------------------------------------------------------------------------#
# Returns Signature specfific features
# For join NMF view_id is the name of one of original matrix to retrieve features from

#' @rdname SignatureSpecificFeatures-methods
#' @aliases SignatureSpecificFeatures,ANY-method
#' @export
#'
#' @examples
#' \dontrun{
#' # For ButchR_joinNMF objects:
#' SignatureSpecificFeatures(jnmf_exp)
#' SignatureSpecificFeatures(jnmf_exp)
#' lapply(SignatureSpecificFeatures(jnmf_exp), function(view){
#'   sapply(view, function(x) sapply(x, length))
#' })
#' lapply(SignatureSpecificFeatures(jnmf_exp, k = 3), function(view){
#'   sapply(view, length)
#' })
#' SignatureSpecificFeatures(jnmf_exp, k = 3, return_all_features = TRUE)
#' SignatureSpecificFeatures(jnmf_exp, k = 3,
#'                           return_all_features = TRUE,
#'                           view_id = "atac")
#' SignatureSpecificFeatures(jnmf_exp,
#'                           return_all_features = TRUE,
#'                           view_id = "atac")
#' }
setMethod("SignatureSpecificFeatures",
          "ButchR_joinNMF",
          function(x, k = NULL, return_all_features = FALSE, view_id = NULL, ...){
            # if no feature extraction was performed
            if (length(x@SignFeatures) == 0) {
              stop("No feature extraction has been performed, please run: \n",
                   "`compute_SignatureFeatures`")
            }

            # Check if view id is indeed one of the views
            if (is.null(view_id)) {
              view_id <- stats::setNames(names(x@SignFeatures), names(x@SignFeatures))
            } else if (all(view_id %in% names(x@SignFeatures))) {
              view_id <- stats::setNames(view_id, view_id)
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




#' @rdname compute_SignatureFeatures-methods
#' @aliases compute_SignatureFeatures,ANY-method
#'
#' @param x An object of class ButchR_joinNMF.
#' @export
#'
#' @examples
#' \dontrun{
#' norm_mat_list <- list(a = matrix(abs(rnorm(1000)), ncol = 10),
#'                       b = matrix(abs(rnorm(1000)), ncol = 10))
#' inmf_exp <- run_joinNMF_tensor(norm_mat_list,
#'                             ranks = 2:5,
#'                             n_initializations     = 10,
#'                             iterations            = 10^4,
#'                             convergence_threshold = 40,
#'                             extract_features = FALSE)
#' inmf_exp <- compute_SignatureFeatures(inmf_exp)
#' SignatureSpecificFeatures(inmf_exp, k = 3)
#' SignatureSpecificFeatures(inmf_exp, k = 3, return_all_features = TRUE)
#' }
setMethod("compute_SignatureFeatures",
          "ButchR_joinNMF",
          function(x){
            k <- x@OptKStats$k
            if (length(k) == 1 & 2 %in% k) {
              stop("Signature Specific Features extraction is not supported for K = 2")
            }
            viewsIDs <- x@input_data$dim$view_ids
            view_specific_WMatrix_list <- x@WMatrix_vs
            #------------------------------------------------------------------#
            #             Compute signatures specific features                 #
            #------------------------------------------------------------------#
            SignFeatures <- lapply(stats::setNames(viewsIDs, viewsIDs), function(viewsID){
              SignFeatures_eval <- lapply(view_specific_WMatrix_list, function(k_eval){
                W <- k_eval[[viewsID]]
                if (ncol(W) == 2) {
                  return(NULL)
                } else {
                  rownames(W) <- x@input_data$rownames[[viewsID]]
                  return(WcomputeFeatureStats(W))
                }
              })
              data.frame(do.call(cbind, SignFeatures_eval),
                         stringsAsFactors = FALSE)
            })


            #------------------------------------------------------------------#
            #                  Return ButchR_joinNMF object                    #
            #------------------------------------------------------------------#
            x@SignFeatures <- SignFeatures
            return(x)
          }
)





