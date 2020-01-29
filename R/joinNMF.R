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
#setGeneric("HMatrix", function(x, k = NULL, ...) standardGeneric("HMatrix"))

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
#setGeneric("WMatrix", function(x, k = NULL, ...) standardGeneric("WMatrix"))

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
#'
#'
#'
#' # run_join_NMF_tensor <- function (matrix_list,
#' #                                  k_min = 2,
#' #                                  k_max = 2,
#' #                                  outer_iter = 10,
#' #                                  inner_iter = 10^4,
#' #                                  conver_stop_threshold = 40,
#' #                                  Sp = 0){
#' #   # Convert params to integer
#' #   nmf_params <- lapply(list(k_min = k_min,
#' #                             k_max = k_max,
#' #                             outer_iter = outer_iter,
#' #                             inner_iter = inner_iter,
#' #                             conver_stop_threshold = conver_stop_threshold),
#' #                        as.integer)
#' #   viewsIDs <- setNames(names(matrix_list), names(matrix_list))
#' #   # Source NMF tensorflow python script
#' #   source_python(file.path(system.file(package = "Bratwurst"), "python/tensor_jNMF.py"))
#' #
#' #   cat("Running join NMF for views: ", paste(names(matrix_list), collapse = ","), "\n")
#' #   #----------------------------------------------------------------------------#
#' #   #     Run integrative NMF - returns list with all ks and all iterations      #
#' #   #----------------------------------------------------------------------------#
#' #   complete_eval <- lapply(k_min:k_max, function(k) {
#' #     k <- as.integer(k)
#' #
#' #     print(Sys.time())
#' #     cat("Factorization rank: ", k, "\n")
#' #     k_eval <- lapply(1:outer_iter, function(i) {
#' #       if (i%%10 == 0) cat("\tIteration: ", i, "\n")
#' #
#' #       jnmf_eval <- jNMF_tensor_py(unname(matrix_list),
#' #                                   rank           = k,
#' #                                   iterations     = nmf_params$inner_iter,
#' #                                   Sp             = Sp,
#' #                                   stop_threshold = nmf_params$conver_stop_threshold)
#' #
#' #       names(jnmf_eval) <- c("Ws", "sharedH", "iterations", "Frob_error")
#' #       names(jnmf_eval$Ws)         <- names(matrix_list)
#' #       names(jnmf_eval$Frob_error) <- names(matrix_list)
#' #       return(jnmf_eval)
#' #     })
#' #     names(k_eval) <- paste0("iter", 1:outer_iter)
#' #
#' #     iters <- paste(sapply(k_eval, function(x) {x$iterations}), collapse = ",")
#' #     print(paste("NMF converged after ", iters, "iterations"))
#' #
#' #     return(k_eval)
#' #   })
#' #
#' #   # Build NMF object
#' #   names(complete_eval) <- paste0("k", k_min:k_max)
#' #
#' #   #----------------------------------------------------------------------------#
#' #   #      Build old NMF experiment objects to return factorization metrics      #
#' #   #----------------------------------------------------------------------------#
#' #   nmfExp_list <- lapply(viewsIDs, function(viewID){
#' #     nmf.exp <- nmfExperimentFromMatrix(matrix_list[[viewID]])
#' #
#' #
#' #     dec.matrix <- lapply(1:(k_max-k_min+1), function(k) {
#' #       k <- as.integer(k)
#' #       k.matrix <- lapply(1:outer_iter, function(i) {
#' #         list(W = complete_eval[[k]][[i]][["Ws"]][[viewID]],
#' #              H = complete_eval[[k]][[i]][["sharedH"]],
#' #              Frob.error = complete_eval[[k]][[i]][["Frob_error"]][[viewID]])
#' #       })
#' #       names(k.matrix) <- 1:outer_iter
#' #       return(k.matrix)
#' #     })
#' #     # Build NMF object
#' #     names(dec.matrix) <- k_min:k_max
#' #     frob.errors <- DataFrame(getFrobError(dec.matrix))
#' #     colnames(frob.errors) <- as.character(k_min:k_max)
#' #     nmf.exp <- setFrobError(nmf.exp, frob.errors)
#' #     nmf.exp <- setHMatrixList(nmf.exp, getHMatrixList(dec.matrix))
#' #     nmf.exp <- setWMatrixList(nmf.exp, getWMatrixList(dec.matrix))
#' #     nmf.exp <- computeFrobErrorStats(nmf.exp)
#' #     nmf.exp <- computeSilhoutteWidth(nmf.exp)
#' #     nmf.exp <- computeCopheneticCoeff(nmf.exp)
#' #     nmf.exp <- computeAmariDistances(nmf.exp)
#' #     return(nmf.exp)
#' #   })
#' #   #----------------------------------------------------------------------------#
#' #   #    Select which outer iteration was the best, based on the Frob error      #
#' #   #----------------------------------------------------------------------------#
#' #   FrobError_list <- lapply(nmfExp_list, function(nmf.exp){
#' #     as.matrix(nmf.exp@FrobError)
#' #   })
#' #   best_factorization_idx <- apply(Reduce("+", FrobError_list), 2, which.min)
#' #   names(best_factorization_idx) <- paste0("k", names(best_factorization_idx))
#' #
#' #   #----------------------------------------------------------------------------#
#' #   #        Organize objects to save to an integrative_NMF class4 object        #
#' #   #----------------------------------------------------------------------------#
#' #   # Shared H matrix list
#' #   shared_HMatrix_list <- lapply(complete_eval, function(k_eval){
#' #     lapply(k_eval, function(jnmf_eval){
#' #       jnmf_eval$sharedH
#' #     })
#' #   })
#' #
#' #   # View specific W matrix list
#' #   view_specific_WMatrix_list <- lapply(complete_eval, function(k_eval){
#' #     lapply(k_eval, function(jnmf_eval){
#' #       jnmf_eval$Ws
#' #     })
#' #   })
#' #
#' #   #----------------------------------------------------------------------------#
#' #   #                       Return integrative_NMF class4 object                 #
#' #   #----------------------------------------------------------------------------#
#' #
#' #   join_NMF(shared_HMatrix_list        = shared_HMatrix_list,
#' #            view_specific_WMatrix_list = view_specific_WMatrix_list,
#' #            best_factorization_idx     = best_factorization_idx,
#' #            view_specific_NMFexp_list  = nmfExp_list)
#' # }
#' #
#'
