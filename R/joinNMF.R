#------------------------------------------------------------------------------#
#                  Join NMF tensorflow Wrapper - CLASS - FUNCTION              #
#------------------------------------------------------------------------------#

#' Join NMF Class
#'
#' @slot shared_HMatrix_list list.
#' @slot view_specific_WMatrix_list list.
#' @slot best_factorization_idx numeric
#' @slot view_specific_NMFexp_list list
#'
#' @return
#'
#' @export
#'
#' @examples
joint_NMF <- setClass(
  Class = "joint_NMF",
  slots = list(shared_HMatrix_list        = "list",
               view_specific_WMatrix_list = "list",
               #frobenius_error            = "data.frame",
               #optimal_K_stats            = "data.frame",
               best_factorization_idx     = "numeric",
               view_specific_NMFexp_list  = "list" )
)

#' Computes integrative NMF on tensorflow using the reticulate framework
#'
#' @param matrix_list
#' @param k_min
#' @param k_max
#' @param outer_iter
#' @param inner_iter
#' @param conver_stop_threshold
#' @param lambda
#'
#' @return
#'
#' @import SummarizedExperiment
#' @importFrom data.table fread
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
run_join_NMF_tensor <- function (matrix_list,
                                 k_min = 2,
                                 k_max = 2,
                                 outer_iter = 10,
                                 inner_iter = 10^4,
                                 conver_stop_threshold = 40,
                                 Sp = 0){
  # Convert params to integer
  nmf_params <- lapply(list(k_min = k_min,
                            k_max = k_max,
                            outer_iter = outer_iter,
                            inner_iter = inner_iter,
                            conver_stop_threshold = conver_stop_threshold),
                       as.integer)
  viewsIDs <- setNames(names(matrix_list), names(matrix_list))
  # Source NMF tensorflow python script
  source_python(file.path(system.file(package = "Bratwurst"), "python/tensor_jNMF.py"))

  cat("Running join NMF for views: ", paste(names(matrix_list), collapse = ","), "\n")
  #----------------------------------------------------------------------------#
  #     Run integrative NMF - returns list with all ks and all iterations      #
  #----------------------------------------------------------------------------#
  complete_eval <- lapply(k_min:k_max, function(k) {
    k <- as.integer(k)

    print(Sys.time())
    cat("Factorization rank: ", k, "\n")
    k_eval <- lapply(1:outer_iter, function(i) {
      if (i%%10 == 0) cat("\tIteration: ", i, "\n")

      jnmf_eval <- jNMF_tensor_py(unname(matrix_list),
                                  rank           = k,
                                  iterations     = nmf_params$inner_iter,
                                  Sp             = Sp,
                                  stop_threshold = nmf_params$conver_stop_threshold)

      names(jnmf_eval) <- c("Ws", "sharedH", "iterations", "Frob_error")
      names(jnmf_eval$Ws)         <- names(matrix_list)
      names(jnmf_eval$Frob_error) <- names(matrix_list)
      return(jnmf_eval)
    })
    names(k_eval) <- paste0("iter", 1:outer_iter)

    iters <- paste(sapply(k_eval, function(x) {x$iterations}), collapse = ",")
    print(paste("NMF converged after ", iters, "iterations"))

    return(k_eval)
  })

  # Build NMF object
  names(complete_eval) <- paste0("k", k_min:k_max)

  #----------------------------------------------------------------------------#
  #      Build old NMF experiment objects to return factorization metrics      #
  #----------------------------------------------------------------------------#
  nmfExp_list <- lapply(viewsIDs, function(viewID){
    nmf.exp <- nmfExperimentFromMatrix(matrix_list[[viewID]])


    dec.matrix <- lapply(1:(k_max-k_min+1), function(k) {
      k <- as.integer(k)
      k.matrix <- lapply(1:outer_iter, function(i) {
        list(W = complete_eval[[k]][[i]][["Ws"]][[viewID]],
             H = complete_eval[[k]][[i]][["sharedH"]],
             Frob.error = complete_eval[[k]][[i]][["Frob_error"]][[viewID]])
      })
      names(k.matrix) <- 1:outer_iter
      return(k.matrix)
    })
    # Build NMF object
    names(dec.matrix) <- k_min:k_max
    frob.errors <- DataFrame(getFrobError(dec.matrix))
    colnames(frob.errors) <- as.character(k_min:k_max)
    nmf.exp <- setFrobError(nmf.exp, frob.errors)
    nmf.exp <- setHMatrixList(nmf.exp, getHMatrixList(dec.matrix))
    nmf.exp <- setWMatrixList(nmf.exp, getWMatrixList(dec.matrix))
    nmf.exp <- computeFrobErrorStats(nmf.exp)
    nmf.exp <- computeSilhoutteWidth(nmf.exp)
    nmf.exp <- computeCopheneticCoeff(nmf.exp)
    nmf.exp <- computeAmariDistances(nmf.exp)
    return(nmf.exp)
  })
  #----------------------------------------------------------------------------#
  #    Select which outer iteration was the best, based on the Frob error      #
  #----------------------------------------------------------------------------#
  FrobError_list <- lapply(nmfExp_list, function(nmf.exp){
    as.matrix(nmf.exp@FrobError)
  })
  best_factorization_idx <- apply(Reduce("+", FrobError_list), 2, which.min)
  names(best_factorization_idx) <- paste0("k", names(best_factorization_idx))

  #----------------------------------------------------------------------------#
  #        Organize objects to save to an integrative_NMF class4 object        #
  #----------------------------------------------------------------------------#
  # Shared H matrix list
  shared_HMatrix_list <- lapply(complete_eval, function(k_eval){
    lapply(k_eval, function(jnmf_eval){
      jnmf_eval$sharedH
    })
  })

  # View specific W matrix list
  view_specific_WMatrix_list <- lapply(complete_eval, function(k_eval){
    lapply(k_eval, function(jnmf_eval){
      jnmf_eval$Ws
    })
  })

  #----------------------------------------------------------------------------#
  #                       Return integrative_NMF class4 object                 #
  #----------------------------------------------------------------------------#

  integrative_NMF(shared_HMatrix_list        = shared_HMatrix_list,
                  view_specific_WMatrix_list = view_specific_WMatrix_list,
                  best_factorization_idx     = best_factorization_idx,
                  view_specific_NMFexp_list  = nmfExp_list)
}


