#' @include generics.R
NULL

#' Integrative NMF
#'
#' Computes integrative NMF on tensorflow using the reticulate framework,
#' uses a list of non-negative matrices as input
#'
#' @param matrix_list List of non-negative matrices
#' @param n_initializations Number of initializations to evaluate
#' @param ranks numeric vector with ranks to factorize
#' @param iterations Maximum number of iterations to run for every initialization
#' @param convergence_threshold The factorization stops,
#' if the convergence test is constant for this number of iterations
#'
#' @return A integrative_NMF object,
#' containing a integrative H matrix and one W matrix for each input matrix
#'
#' @export
#'
#' @examples
#' inmf_exp <- run_iNMF_tensor(list(a = matrix(1:1000, ncol = 10),
#'                                     b = matrix(1:1000, ncol = 10)),
#'                                ranks = 2:5,
#'                                n_initializations     = 10,
#'                                iterations            = 10^4,
#'                                convergence_threshold = 40)
#' inmf_exp
run_iNMF_tensor <- function (matrix_list,
                             ranks                 = 2,
                             n_initializations     = 10,
                             iterations            = 10^4,
                             convergence_threshold = 40,
                             Sp                    = 0,
                             lamb                  = 10,
                             extract_features      = FALSE){

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
  #    Run integrative NMF - returns list with all ks and all iterations       #
  #----------------------------------------------------------------------------#
  # Source iNMF tensorflow python function
  iNMF_tensor_py <- source_NMFtensor_function("iNMF")

  # Run jNMF
  cat("Running integrative NMF for views: ", paste(names(matrix_list), collapse = ","), "\n")
  complete_eval <- lapply(nmf_params$ranks, function(k) {
    print(Sys.time())
    cat("Factorization rank: ", k, "\n")
    k_eval <- iNMF_tensor_py(matrix_list       = unname(matrix_list),
                             rank              = k,
                             n_initializations = nmf_params$n_initializations,
                             iterations        = nmf_params$iterations,
                             Sp                = Sp,
                             stop_threshold    = nmf_params$convergence_threshold,
                             lamb              = lamb)

    names(k_eval)     <- c("Ws", "sharedH", "Hs", "iterations", "Frob_error", "W_eval")
    names(k_eval$Ws)  <- viewsIDs
    names(k_eval$Hs)  <- viewsIDs
    k_eval$iterations <- unlist(k_eval$iterations)
    k_eval$Frob_error <- unlist(k_eval$Frob_error)

    # Optimal K stats
    k_eval$OptKStats <- try(compute_OptKStats_NMF(k_eval, k), silent = FALSE)
    k_eval$W_eval <- NULL

    print(paste("integrative NMF converged after ",
                paste(k_eval$iterations, collapse = ","),
                "iterations"))
    return(k_eval)
  })
  #----------------------------------------------------------------------------#
  #              Build integrative NMF object slots                            #
  #----------------------------------------------------------------------------#
  # input data info
  input_data <- list(hash       = digest::digest(matrix_list),
                     dim        = data.frame(view_ids = names(matrix_list),
                                       do.call(rbind, lapply(matrix_list,
                                                             dim))),
                     colnames   = colnames(matrix_list[[1]]),
                     rownames   = lapply(matrix_list, rownames),
                     run_params = list(n_initializations = nmf_params$n_initializations,
                                       iterations        = nmf_params$iterations,
                                       stop_threshold    = nmf_params$convergence_threshold,
                                       Sp                = Sp,
                                       lamb              = lamb,
                                       extract_features  = extract_features))




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


  # View specific W matrix list
  view_specific_WMatrix_list <- lapply(complete_eval, function(k_eval){
    k_eval$Ws
  })

  # Shared H matrix list
  shared_HMatrix_list <- lapply(complete_eval, function(k_eval){
    k_eval$sharedH
  })

  # View specific H matrix list
  view_specific_HMatrix_list <- lapply(complete_eval, function(k_eval){
    k_eval$Hs
  })

  #----------------------------------------------------------------------------#
  #                  Compute signatures specific features                      #
  #----------------------------------------------------------------------------#
  if (extract_features) {
    SignFeatures <- lapply(viewsIDs, function(viewsID){

      #print(names(k_eval))
      SignFeatures_eval <- lapply(view_specific_WMatrix_list, function(k_eval){
        W <- k_eval[[viewsID]]
        if (ncol(W) == 2) {
          return(NULL)
        } else {
          rownames(W) <- input_data$rownames[[viewsID]]
          return(WcomputeFeatureStats(W))
        }
      })
      SignFeatures_eval <- data.frame(do.call(cbind, SignFeatures_eval),
                                      stringsAsFactors = FALSE)
      return(SignFeatures_eval)
    })
  } else {
    SignFeatures <- list()
  }

  #----------------------------------------------------------------------------#
  #                       Return integrative_NMF object                        #
  #----------------------------------------------------------------------------#

  integrative_NMF(input_data   = input_data,
                  HMatrix      = shared_HMatrix_list,
                  HMatrix_vs   = view_specific_HMatrix_list,
                  WMatrix_vs   = view_specific_WMatrix_list,
                  FrobError    = frob_errors,
                  OptKStats    = OptKStats,
                  OptK         = OptK,
                  SignFeatures = SignFeatures)
}




#'
#'
#'
#' #' Computes integrative NMF on tensorflow using the reticulate framework
#' #'
#' #' @param matrix_list
#' #' @param k_min
#' #' @param k_max
#' #' @param outer_iter
#' #' @param inner_iter
#' #' @param conver_stop_threshold
#' #' @param lambda
#' #'
#' #' @return
#' #'
#' #' @import SummarizedExperiment
#' #' @importFrom data.table fread
#' #' @importFrom S4Vectors DataFrame
#' #' @export
#' #'
#' #' @examples
#' run_integrative_NMF_tensor <- function (matrix_list,
#'                                         k_min = 2,
#'                                         k_max = 2,
#'                                         outer_iter = 10,
#'                                         inner_iter = 10^4,
#'                                         conver_stop_threshold = 40,
#'                                         lambda = 1,
#'                                         Sp = 0){
#'   # Convert params to integer
#'   nmf_params <- lapply(list(k_min = k_min,
#'                             k_max = k_max,
#'                             outer_iter = outer_iter,
#'                             inner_iter = inner_iter,
#'                             conver_stop_threshold = conver_stop_threshold,
#'                             lambda = lambda),
#'                        as.integer)
#'   viewsIDs <- setNames(names(matrix_list), names(matrix_list))
#'   # Source NMF tensorflow python script
#'   source_python(file.path(system.file(package = "Bratwurst"), "python/intnmf_tensor.py"))
#'
#'   cat("Running integrative NMF for views: ", paste(names(matrix_list), collapse = ","), "\n")
#'   #============================================================================#
#'   #     Run integrative NMF - returns list with all ks and all iterations      #
#'   #============================================================================#
#'   complete_eval <- lapply(k_min:k_max, function(k) {
#'     k <- as.integer(k)
#'
#'     print(Sys.time())
#'     cat("Factorization rank: ", k, "\n")
#'     k_eval <- lapply(1:outer_iter, function(i) {
#'       if (i%%10 == 0) cat("\tIteration: ", i, "\n")
#'
#'       inmf_eval <- iNMF(unname(matrix_list),
#'                         rank           = k,
#'                         iterations     = nmf_params$inner_iter,
#'                         L              = nmf_params$lambda,
#'                         Sp             = Sp,
#'                         stop_threshold = nmf_params$conver_stop_threshold)
#'
#'       names(inmf_eval) <- c("Ws", "sharedH", "viewHs", "iterations", "Frob_error")
#'       names(inmf_eval$Ws)         <- names(matrix_list)
#'       names(inmf_eval$viewHs)     <- names(matrix_list)
#'       names(inmf_eval$Frob_error) <- names(matrix_list)
#'       return(inmf_eval)
#'     })
#'     names(k_eval) <- paste0("iter", 1:outer_iter)
#'
#'     iters <- paste(sapply(k_eval, function(x) {x$iterations}), collapse = ",")
#'     print(paste("NMF converged after ", iters, "iterations"))
#'
#'     return(k_eval)
#'   })
#'
#'   # Build NMF object
#'   names(complete_eval) <- paste0("k", k_min:k_max)
#'
#'   #============================================================================#
#'   #      Build old NMF experiment objects to return factorization metrics      #
#'   #============================================================================#
#'   nmfExp_list <- lapply(viewsIDs, function(viewID){
#'     nmf.exp <- nmfExperimentFromMatrix(matrix_list[[viewID]])
#'
#'
#'     dec.matrix <- lapply(1:(k_max-k_min+1), function(k) {
#'       k <- as.integer(k)
#'       k.matrix <- lapply(1:outer_iter, function(i) {
#'         list(W = complete_eval[[k]][[i]][["Ws"]][[viewID]],
#'              H = complete_eval[[k]][[i]][["sharedH"]] + complete_eval[[k]][[i]][["viewHs"]][[viewID]],
#'              Frob.error = complete_eval[[k]][[i]][["Frob_error"]][[viewID]])
#'       })
#'       names(k.matrix) <- 1:outer_iter
#'       return(k.matrix)
#'     })
#'     # Build NMF object
#'     names(dec.matrix) <- k_min:k_max
#'     frob.errors <- DataFrame(getFrobError(dec.matrix))
#'     colnames(frob.errors) <- as.character(k_min:k_max)
#'     nmf.exp <- setFrobError(nmf.exp, frob.errors)
#'     nmf.exp <- setHMatrixList(nmf.exp, getHMatrixList(dec.matrix))
#'     nmf.exp <- setWMatrixList(nmf.exp, getWMatrixList(dec.matrix))
#'     nmf.exp <- computeFrobErrorStats(nmf.exp)
#'     nmf.exp <- computeSilhoutteWidth(nmf.exp)
#'     nmf.exp <- computeCopheneticCoeff(nmf.exp)
#'     nmf.exp <- computeAmariDistances(nmf.exp)
#'     return(nmf.exp)
#'   })
#'   #============================================================================#
#'   #    Select which outer iteration was the best, based on the Frob error      #
#'   #============================================================================#
#'   FrobError_list <- lapply(nmfExp_list, function(nmf.exp){
#'     as.matrix(nmf.exp@FrobError)
#'   })
#'   best_factorization_idx <- apply(Reduce("+", FrobError_list), 2, which.min)
#'   names(best_factorization_idx) <- paste0("k", names(best_factorization_idx))
#'
#'   #============================================================================#
#'   #        Organize objects to save to an integrative_NMF class4 object        #
#'   #============================================================================#
#'   # Shared H matrix list
#'   shared_HMatrix_list <- lapply(complete_eval, function(k_eval){
#'     lapply(k_eval, function(inmf_eval){
#'       inmf_eval$sharedH
#'     })
#'   })
#'   # View specific H matrix list
#'   view_specific_HMatrix_list <- lapply(complete_eval, function(k_eval){
#'     lapply(k_eval, function(inmf_eval){
#'       inmf_eval$viewHs
#'     })
#'   })
#'   # View specific W matrix list
#'   view_specific_WMatrix_list <- lapply(complete_eval, function(k_eval){
#'     lapply(k_eval, function(inmf_eval){
#'       inmf_eval$Ws
#'     })
#'   })
#'
#'   #============================================================================#
#'   #                       Return integrative_NMF class4 object                 #
#'   #============================================================================#
#'
#'   integrative_NMF(shared_HMatrix_list        = shared_HMatrix_list,
#'                   view_specific_HMatrix_list = view_specific_HMatrix_list,
#'                   view_specific_WMatrix_list = view_specific_WMatrix_list,
#'                   best_factorization_idx     = best_factorization_idx,
#'                   view_specific_NMFexp_list  = nmfExp_list)
#' }