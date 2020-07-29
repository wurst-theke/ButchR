#' @include generics.R
NULL

#' Integrative NMF
#'
#' Computes integrative NMF on tensorflow using the reticulate framework,
#' uses a list of non-negative matrices as input
#'
#' @param matrix_list List of non-negative matrices
#' @param ranks numeric vector with ranks to factorize
#' @param n_initializations Number of initializations to evaluate
#' @param iterations Maximum number of iterations to run for every initialization
#' @param convergence_threshold The factorization stops,
#' if the convergence test is constant for this number of iterations
#' @param Sp Sparsity,
#' @param lamb Free parameter lambda bigger values shift the decomposition
#' towards finding the common effect between niews.
#' @param extract_features if TRUE performs feature extraction for all
#' factorization ranks > 2.
#'
#' @return An object of class ButchR_integrativeNMF.
#' containing a integrative H matrix and one W matrix for each input matrix
#' @import reticulate
#' @export
#'
#' @examples
#' \dontrun{
#' inmf_exp <- run_iNMF_tensor(list(a = matrix(1:1000, ncol = 10),
#'                                     b = matrix(1:1000, ncol = 10)),
#'                                ranks = 2:5,
#'                                n_initializations     = 10,
#'                                iterations            = 10^4,
#'                                convergence_threshold = 40)
#' inmf_exp
#' }
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
  if (!is.list(matrix_list) ) {
    stop("\nmatrix_list should be a list of Non-negative matrices\n")
  }
  matrix_list <- lapply(matrix_list, val_nonnegative_matrix)
  if (!length(unique(sapply(matrix_list, ncol))) == 1 ) {
    stop("\nThe number of columns should be identical for all matrices\n")
  }
  if (!all(sapply(lapply(matrix_list, colnames),
                  identical, colnames(matrix_list[[1]])))) {
    stop("\nColumn names should be identical between matrices\n")
  }
  if (is.null(names(matrix_list))) {
    names(matrix_list) <- paste0("view", 1:length(matrix_list))
    warning("Input matrix list do not have names, assigning ids:\n",
            paste0(names(matrix_list), collapse = ","), "\n")
  }
  val_ranks_torun(ranks, ncol(matrix_list[[1]]))
  val_single_integer(n_initializations, "n_initializations")
  val_single_integer(iterations, "iterations")
  val_single_integer(convergence_threshold, "convergence_threshold")
  val_single_numeric(Sp, "Sp")
  val_single_numeric(lamb, "lamb")
  val_single_boolean(extract_features, "extract_features")

  # Convert params to integer
  nmf_params <- lapply(list(ranks                 = ranks,
                            n_initializations     = n_initializations,
                            iterations            = iterations,
                            convergence_threshold = convergence_threshold),
                       as.integer)
  names(nmf_params$ranks) <- paste0("k", nmf_params$ranks)
  viewsIDs <- stats::setNames(names(matrix_list), names(matrix_list))
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
  if (length(ranks) == 1 & 2 %in% ranks) {
    SignFeatures <- list()
  }
  #----------------------------------------------------------------------------#
  #                   Return ButchR_integrativeNMF object                      #
  #----------------------------------------------------------------------------#
  ButchR_integrativeNMF(input_data   = input_data,
                        HMatrix      = shared_HMatrix_list,
                        HMatrix_vs   = view_specific_HMatrix_list,
                        WMatrix_vs   = view_specific_WMatrix_list,
                        FrobError    = frob_errors,
                        OptKStats    = OptKStats,
                        OptK         = OptK,
                        SignFeatures = SignFeatures)
}







#' Integrative NMF Tune Lambda
#'
#' Based on Yang and Michailidis, 2016,
#' to tune the value of the parameter lambda for the integrative NMF (iNMF),
#' the objectives values of join NMF (jNMF) are compared to
#' single view NMFs (sNMF), the principle is that join NMF
#' represents complete homogeneity and single view NMF
#' represents complete heterogeneity.
#' To avoid overfitting the best lambda can be selected by plotting
#' the difference of the unsquared residual quantities of
#' jNMF and iNMF (Ri - Rj) over multiple values of lambda,
#' and compare it to the difference of the unsquared residual quantities of
#' sNMF and jNMF c*(Rj - Rs).
#' The optimal lambda usually is the first lambda in which
#' (Ri - Rj) < c*(Rj - Rs).
#' Where c is a constant >= 2.
#'
#' @param matrix_list List of non-negative matrices.
#' @param lambdas a sequence of lambdas to test.
#' @param Output_type Type of desired output, could be:
#' \itemize{
#' \item residuals - Residual quantities data.frame.
#' \item iNMF - Returns the integrative NMF object for the optimal estimated lambda.
#' \item all_iNMF - Returns a list with an integrative NMF object for every  lambda.
#' \item plot - Residual quantities plot,
#' (Ri - Rj) and 2(Rj - Rs) vs lamda range.
#' \item all - list containing all computed NMF objects,
#' Residual quantities data.frame and Residual quantities plot.
#' }
#' @param thr_cons numeric value, Threshold constant c, in  c*(Rj - Rs).
#' @param rank numeric vector with rank to factorize.
#' @param n_initializations Number of initializations to evaluate.
#' @param iterations Maximum number of iterations to run for every initialization.
#' @param convergence_threshold The factorization stops,
#' if the convergence test is constant for this number of iterations.
#' @param Sp Sparcity constrain, values > 0 force sparcity in the H matrix.
#' @param show_plot \code{TRUE} print plot;
#' \code{FALSE} plot is only returned if it is selected in Output_type.
#' @param extract_features \code{TRUE} Extract signature specific features.
#'
#'
#' @return There are five different types of possible outputs, depending on the
#' selected option on `Output_type`. The default output is a data.frame with
#' the residual resulting from tuning lambdas across multiple lambdas.
#' @import ggplot2 dplyr
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare the residual across multiple lambdas:
#' iNMF_lambda_tuning(matrix_list           = norm_mat_list,
#'                    lambdas               = seq(0, 1, 0.1),
#'                    Output_type           = "residuals",
#'                    rank                  = 9,
#'                    n_initializations     = 5,
#'                    iterations            = 10^4,
#'                    convergence_threshold = 40,
#'                    Sp                    = 0,
#'                    extract_features      = FALSE)
#' # Retrieve all the objects and extract the iNMF for the best lambda:
#' inmf_tune <- iNMF_lambda_tuning(matrix_list           = norm_mat_list,
#'                                 lambdas               = seq(0, 1, 0.1),
#'                                 thr_cons              = 4,
#'                                 Output_type           = "all",
#'                                 rank                  = 9,
#'                                 n_initializations     = 5,
#'                                 iterations            = 10^4,
#'                                 convergence_threshold = 40)
#' min(inmf_tune$residuals$lambda[inmf_tune$residuals$best_lambda])
#' inmf_tune$iNMF$lambda_0.2
#' }
iNMF_lambda_tuning <- function (matrix_list,
                                lambdas               = seq(0, 1, 0.1),
                                Output_type           = "residuals",
                                thr_cons              = 3,
                                rank                  = 2,
                                n_initializations     = 10,
                                iterations            = 10^4,
                                convergence_threshold = 40,
                                Sp                    = 0,
                                show_plot             = TRUE,
                                extract_features      = FALSE){

  #----------------------------------------------------------------------------#
  #                            Setup and data check                            #
  #----------------------------------------------------------------------------#
  # Check  data
  if (!is.list(matrix_list) ) {
    stop("\nmatrix_list should be a list of Non-negative matrices\n")
  }
  matrix_list <- lapply(matrix_list, val_nonnegative_matrix)
  if (!length(unique(sapply(matrix_list, ncol))) == 1 ) {
    stop("\nThe number of columns should be identical for all matrices\n")
  }
  if (!all(sapply(lapply(matrix_list, colnames),
                  identical, colnames(matrix_list[[1]])))) {
    stop("\nColumn names should be identical between matrices\n")
  }
  if (is.null(names(matrix_list))) {
    names(matrix_list) <- paste0("view", 1:length(matrix_list))
    warning("Input matrix list do not have names, assigning ids:\n",
            paste0(names(matrix_list), collapse = ","), "\n")
  }
  val_positive_numeric(lambdas, "lambdas")
  if (!is.character(Output_type)) {
    stop("\nOutput_type invalid value, output types supported are only:\n
         'residuals', 'iNMF', 'all_iNMF', 'plot' or 'all' \n")
  }
  if (length(Output_type) != 1  |
      !any(Output_type %in% c("residuals", "iNMF", "all_iNMF", "plot", "all"))) {
    stop("\nOutput_type invalid value, output types supported are only:\n
         'residuals', 'iNMF', 'all_iNMF', 'plot' or 'all' \n")
  }
  # if (!any(Output_type %in% c("residuals", "iNMF", "all_iNMF", "plot", "all")) ) {
  #   stop("\nOutput_type invalid value, output types supported are only:\n
  #        'residuals', 'iNMF', 'all_iNMF', 'plot' or 'all' \n")
  # }
  if (!is.logical(show_plot) & !length(show_plot) == 1) {
    stop("\nshow_plot invalid value, select TRUE or FALSE \n")
  }
  if (!length(rank) == 1) {
    stop("\nrank has to be a unique integer number \n")
  }






  # val_ranks_torun(ranks, ncol(matrix_list[[1]]))
  # val_single_integer(n_initializations, "n_initializations")
  # val_single_integer(iterations, "iterations")
  # val_single_integer(convergence_threshold, "convergence_threshold")
  # val_single_numeric(Sp, "Sp")
  # val_single_numeric(lamb, "lamb")
  # val_single_boolean(extract_features, "extract_features")




  # Convert params to integer
  nmf_params <- lapply(list(rank                  = rank,
                            n_initializations     = n_initializations,
                            iterations            = iterations,
                            convergence_threshold = convergence_threshold),
                       as.integer)
  names(nmf_params$rank) <- paste0("k", nmf_params$rank)
  viewsIDs <- stats::setNames(names(matrix_list), names(matrix_list))
  #----------------------------------------------------------------------------#
  #                                    Join NMF                                #
  #----------------------------------------------------------------------------#
  jnmf_obj <- run_joinNMF_tensor(matrix_list           = matrix_list,
                                 ranks                 = nmf_params$rank,
                                 n_initializations     = nmf_params$n_initializations,
                                 iterations            = nmf_params$iterations,
                                 convergence_threshold = nmf_params$convergence_threshold,
                                 Sp                    = Sp,
                                 extract_features      = extract_features)
  # Estimate Frob. norm
  norm_jnmf <- lapply(viewsIDs, function(viewID){
    #xn <- norm(matrix_list[[viewID]], type = "F")
    x <- matrix_list[[viewID]]
    w <- WMatrix(jnmf_obj, k = nmf_params$rank, view_id = viewID)
    h <- HMatrix(jnmf_obj, k = nmf_params$rank)
    #fn <- norm((matrix_list[[viewID]] - (w %*% h)), type = "F")
    norm((x - (w %*% h)), type = "F")
  })
  norm_jnmf <- sum(unlist(norm_jnmf))

  #----------------------------------------------------------------------------#
  #                                Single view NMF                             #
  #----------------------------------------------------------------------------#
  snmf_obj_list <- lapply(viewsIDs, function(viewID){
    cat("Running single view NMF for : ", viewID, "\n")
    run_NMF_tensor(X = matrix_list[[viewID]],
                      ranks                 = nmf_params$rank,
                      method = "NMF",
                      n_initializations     = nmf_params$n_initializations,
                      iterations            = nmf_params$iterations,
                      convergence_threshold = nmf_params$convergence_threshold,
                      extract_features      = extract_features)
  })
  # Estimate Frob. norm
  norm_snmf <- lapply(viewsIDs, function(viewID){
    x <- matrix_list[[viewID]]
    w <- WMatrix(snmf_obj_list[[viewID]], k = nmf_params$rank)
    h <- HMatrix(snmf_obj_list[[viewID]], k = nmf_params$rank)
    norm((x - (w %*% h)), type = "F")
  })
  norm_snmf <- sum(unlist(norm_snmf))


  #thres <- norm_jnmf + c*(norm_jnmf - norm_snmf)
  #----------------------------------------------------------------------------#
  #                              Integrative NMF                               #
  #----------------------------------------------------------------------------#
  names(lambdas) <- paste0("lambda_", lambdas)
  inmf_obj_list <- lapply(lambdas, function(lamb){
    cat("Testing lamda : ", lamb, "\n")
    inmf_obj <- run_iNMF_tensor(matrix_list           = matrix_list,
                                ranks                 = nmf_params$rank,
                                n_initializations     = nmf_params$n_initializations,
                                iterations            = nmf_params$iterations,
                                convergence_threshold = nmf_params$convergence_threshold,
                                Sp                    = Sp,
                                lamb                  = lamb,
                                extract_features      = extract_features)
  })
  # Estimate Frob. norm
  norm_inmf <- sapply(inmf_obj_list, function(inmf_obj){
    norm_inmf_sv <- lapply(viewsIDs, function(viewID){
      x <- matrix_list[[viewID]]
      w <- WMatrix(inmf_obj, k = nmf_params$rank, view_id = viewID)
      h <- HMatrix(inmf_obj, k = nmf_params$rank, type = "shared")
      norm((x - (w %*% h)), type = "F")
    })
    norm_inmf_sv <- sum(unlist(norm_inmf_sv))
  })

  residuals_df <- data.frame(norm_jnmf = norm_jnmf,
                             norm_snmf = norm_snmf,
                             #thres     = thres,
                             norm_inmf = norm_inmf,
                             lambda    = lambdas,
                             stringsAsFactors = FALSE)
  residuals_df <- residuals_df %>%
    dplyr::mutate(diff_iNMF_jNMF = norm_inmf - norm_jnmf) %>%
    dplyr::mutate(diff_iNMF_sNMF = abs(thr_cons*(norm_jnmf - norm_snmf))) %>%
    #dplyr::mutate(diff_iNMF_sNMF = thr_cons*(norm_jnmf - norm_snmf)) %>%
    dplyr::arrange(-.data$lambda) %>%
    dplyr::mutate(best_lambda = cumsum(.data$diff_iNMF_jNMF > .data$diff_iNMF_sNMF) == 1) %>%
    dplyr::mutate(best_lambda = c(.data$best_lambda[-1], .data$best_lambda[1]))

  if (Output_type %in% c("plot", "all") | show_plot) {
    residuals_gg <- residuals_df %>%
      dplyr::select(.data$lambda, .data$diff_iNMF_jNMF,
                    .data$diff_iNMF_sNMF, .data$best_lambda) %>%
      tidyr::pivot_longer(-c("lambda", "best_lambda"),
                          names_to = "norm",
                          values_to = "unsquared_residual_quantities") %>%
      ggplot(aes(x = .data$lambda, y = .data$unsquared_residual_quantities, color = norm)) +
      geom_vline(data = function(x){x %>% dplyr::filter(.data$best_lambda)},
                 aes(xintercept = .data$lambda)) +
      geom_line() +
      theme_bw()
    if (show_plot & !Output_type == "plot") print(residuals_gg)
  }

  if (Output_type == "plot") {
    return(residuals_gg)
  } else if (Output_type == "residuals") {
    return(residuals_df)
  } else if (Output_type == "iNMF") {
    idx <- which(lambdas == min(residuals_df$lambda[residuals_df$best_lambda]))
    return(inmf_obj_list[[idx]])
  } else if (Output_type == "all_iNMF") {
    return(inmf_obj_list)
  } else if (Output_type == "all") {
    return(list(iNMF      = inmf_obj_list,
                jNMF      = jnmf_obj,
                sNMF      = snmf_obj_list,
                plot      = residuals_gg,
                residuals = residuals_df))
  }
}



