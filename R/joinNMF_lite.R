#' @include generics.R
NULL

#' Join NMF
#'
#' Computes join NMF on tensorflow using the reticulate framework,
#' uses a list of non-negative matrices as input
#'
#' @param matrix_list List of non-negative matrices
#' @param n_initializations Number of initializations to evaluate
#' @param ranks numeric vector with ranks to factorize
#' @param iterations Maximum number of iterations to run for every initialization
#' @param convergence_threshold The factorization stops,
#' if the convergence test is constant for this number of iterations
#'
#' @return A join_NMF object,
#' containing a join H matrix and one W matrix for each input matrix
#'
#' @export
#'
#' @examples
#' jnmf_exp <- run_joinNMF_tensor(list(a = matrix(1:1000, ncol = 10),
#'                                     b = matrix(1:1000, ncol = 10)),
#'                                ranks = 2:5,
#'                                n_initializations     = 10,
#'                                iterations            = 10^4,
#'                                convergence_threshold = 40)
#' jnmf_exp
run_joinNMF_tensor <- function (matrix_list,
                                ranks  = 2,
                                n_initializations     = 10,
                                iterations            = 10^4,
                                convergence_threshold = 40,
                                Sp = 0){

  #----------------------------------------------------------------------------#
  #                            Setup and data check                            #
  #----------------------------------------------------------------------------#
  # Check  data
  if (!is.list(matrix_list) &
      (!sum(sapply(matrix_list, is.matrix)) == length(matrix_list))) {
    stop("\nmatrix_list should be a list of Non-negative matrices\n")
  }
  if (min(sapply(matrix_list, min)) < 0 ) {
    stop("\nNegative values present in input matrix\n
         only non-negative matrices supported\n")
  }
  if (!all(sapply(lapply(matrix_list, colnames),
                  identical, colnames(matrix_list[[1]])))) {
    stop("\nColumn names should be identical between matrices\n")
  }
  if (is.null(names(matrix_list))) {
    names(matrix_list) <- paste0("view", 1:length(matrix_list))
    warning("Input matrix list do not have names, assigning ids:\n", names(matrix_list), "\n")
  }

  # Convert params to integer
  nmf_params <- lapply(list(ranks                 = ranks,
                            n_initializations     = n_initializations,
                            iterations            = iterations,
                            convergence_threshold = convergence_threshold),
                       as.integer)
  names(nmf_params$ranks) <- paste0("k", nmf_params$ranks)
  viewsIDs <- setNames(names(matrix_list), names(matrix_list))
  #----------------------------------------------------------------------------#
  #          Run join NMF - returns list with all ks and all iterations        #
  #----------------------------------------------------------------------------#
  # Source jNMF tensorflow python function
  jNMF_tensor_py <- source_NMFtensor_function("jNMF")

  # Run jNMF
  cat("Running join NMF for views: ", paste(names(matrix_list), collapse = ","), "\n")
  complete_eval <- lapply(nmf_params$ranks, function(k) {
    print(Sys.time())
    cat("Factorization rank: ", k, "\n")
    k_eval <- jNMF_tensor_py(unname(matrix_list),
                             rank              = k,
                             n_initializations = nmf_params$n_initializations,
                             iterations        = nmf_params$iterations,
                             Sp                = Sp,
                             stop_threshold    = nmf_params$convergence_threshold)

    names(k_eval)     <- c("Ws", "sharedH", "iterations", "Frob_error", "W_eval")
    names(k_eval$Ws)  <- viewsIDs
    k_eval$iterations <- unlist(k_eval$iterations)
    k_eval$Frob_error <- unlist(k_eval$Frob_error)

    # Optimal K stats
    k_eval$OptKStats <- try(compute_OptKStats_NMF(k_eval, k), silent = FALSE)
    k_eval$W_eval <- NULL

    print(paste("join NMF converged after ",
                paste(k_eval$iterations, collapse = ","),
                "iterations"))
    return(k_eval)
  })
  #----------------------------------------------------------------------------#
  #                    Build join NMF object slots                             #
  #----------------------------------------------------------------------------#
  # input data info
  input_data <- list(hash = digest::digest(matrix_list),
                     dim  = data.frame(view_ids = names(matrix_list),
                                       do.call(rbind, lapply(matrix_list,
                                                             dim))),
                     colnames = colnames(matrix_list[[1]]),
                     rownames = lapply(matrix_list, rownames))

  # Frob. error data frame
  frob_errors <- as.data.frame(do.call(cbind, lapply(complete_eval, "[[" , "Frob_error")))


  # Optimal K stats
  OptKStats <- lapply(complete_eval, "[[" , "OptKStats")
  if (!any(sapply(OptKStats, inherits, "try-error"))) {
    OptKStats <- as.data.frame(dplyr::bind_rows(OptKStats))

    # Optimal K
    indCopheneticCoeff <- which(local.maxima(OptKStats$copheneticCoeff)) # Max Cophenetic Coeff
    indMeanAmariDist   <- which(local.minima(OptKStats$meanAmariDist))   # Min Amari Dist
    OptK <- OptKStats$k[intersect(indCopheneticCoeff, indMeanAmariDist)]
    if (length(OptK) == 0) {
      #warning("No optimal K could be determined from the Optimal K stat\n")
      cat("No optimal K could be determined from the Optimal K stat\n")
    }

  } else {
    OptKStats <- data.frame()
    OptK <- integer()
    cat("Error found while computing factorization stats\nSkipping Optimal K\n")
  }

  # Shared H matrix list
  shared_HMatrix_list <- lapply(complete_eval, function(k_eval){
    k_eval$sharedH
  })

  # View specific W matrix list
  view_specific_WMatrix_list <- lapply(complete_eval, function(k_eval){
    k_eval$Ws
  })

  #----------------------------------------------------------------------------#
  #                              Return join_NMF object                        #
  #----------------------------------------------------------------------------#

  join_NMF(input_data = input_data,
           HMatrix    = shared_HMatrix_list,
           WMatrix_vs = view_specific_WMatrix_list,
           FrobError  = frob_errors,
           OptKStats  = OptKStats,
           OptK       = OptK)
}



# #
# #
# environment(run_joinNMF_tensor) <- asNamespace("Bratwurst")
# jnmf_exp <- run_joinNMF_tensor(norm_mat_list,
#                                 ranks  = 4:10,
#                                 n_initializations     = 3,
#                                 iterations            = 10^4,
#                                 convergence_threshold = 40,
#                                 Sp = 0)
# gg_plotKStats(jnmf_exp)
# HMatrix(jnmf_exp)

#
# jnmf_exp
# FrobError(jnmf_exp)
# HMatrix(jnmf_exp)
# lapply(WMatrix(jnmf_exp)[[1]], head)
# jnmf_exp@OptKStats
# jnmf_exp@OptK
#
# all.equal(lapply(norm_mat_list, colnames))
# matrix_list
#
# x <- list(a = matrix(1:1000, ncol = 10),
#           b = matrix(1:1000, ncol = 10))
#
# x <- try(aa, silent = TRUE)
# inherits(x, "try-error")
