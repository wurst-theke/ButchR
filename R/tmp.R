#' runNMFtensor_lite_eval <- function (X,
#'                                                       ranks,
#'                                                       n_initializations     = 10,
#'                                                       iterations            = 10^4,
#'                                                       convergence_threshold = 40){
#'
#'   #----------------------------------------------------------------------------#
#'   #                                    Setup                                   #
#'   #----------------------------------------------------------------------------#
#'   # Check  data
#'   if (min(X) < 0 ) {
#'     stop("\nNegative values present in input matrix\n
#'          only non-negative matrices supported\n")
#'   }
#'   # Convert params to integer
#'   nmf_params <- lapply(list(ranks                 = ranks,
#'                             n_initializations     = n_initializations,
#'                             iterations            = iterations,
#'                             convergence_threshold = convergence_threshold),
#'                        as.integer)
#'
#'   #----------------------------------------------------------------------------#
#'   #                Run single view NMF on tensorflow                           #
#'   #----------------------------------------------------------------------------#
#'   # Source NMF tensorflow python script
#'   source_python(file.path(system.file(package = "Bratwurst"), "python/nmf_tensor_lite.py"))
#'
#'   # Run NMF
#'   #X <- input_matrix
#'   names(nmf_params$ranks) <- paste0("k", nmf_params$ranks)
#'   complete_eval <- lapply(nmf_params$ranks, function(k) {
#'     #k <- as.integer(k)
#'
#'     print(Sys.time())
#'     cat("Factorization rank: ", k, "\n")
#'
#'     k_eval <- NMF_tensor_py(matrix            = X,
#'                             rank              = k,
#'                             n_initializations = nmf_params$n_initializations,
#'                             iterations        = nmf_params$iterations,
#'                             stop_threshold    = nmf_params$convergence_threshold)
#'     names(k_eval) <- c("W", "H", "iterations", "Frob_error", "W_eval")
#'     k_eval$iterations <- unlist(k_eval$iterations)
#'     k_eval$Frob_error <- unlist(k_eval$Frob_error)
#'
#'
#'     print(paste("NMF converged after ", paste(k_eval$iterations, collapse = ","), "iterations"))
#'     return(k_eval)
#'   })
#'
#'   #----------------------------------------------------------------------------#
#'   #                        Build NMF object slots                              #
#'   #----------------------------------------------------------------------------#
#'   # input matrix info
#'   input_matrix <- list(hash     = digest::digest(X),
#'                        dim      = dim(X),
#'                        colnames = colnames(X),
#'                        rownames = rownames(X))
#'
#'   # Frob. error data frame
#'   frob_errors <- as.data.frame(do.call(cbind, lapply(complete_eval, "[[" , "Frob_error")))
#'
#'   # Optimal K stats
#'   OptKStats <- compute_OptKStats_NMF(complete_eval)
#'
#'   # Optimal K
#'   indCopheneticCoeff <- which(local.maxima(OptKStats$copheneticCoeff)) # Max Cophenetic Coeff
#'   indMeanAmariDist   <- which(local.minima(OptKStats$meanAmariDist))   # Min Amari Dist
#'   print(indCopheneticCoeff)
#'   print(indMeanAmariDist)
#'   OptK <- OptKStats$k[intersect(indCopheneticCoeff, indMeanAmariDist)]
#'   if (length(OptK) == 0) {
#'     warning("No optimal K could be determined from the Optimal K stat\n")
#'   }
#'
#'   #print(frob_errors)
#'   #sapply(complete_eval, function(x) print(x$Frob_error))
#'   #print(sapply(complete_eval, function(x) x$Frob_error))
#'   # Factorization_ranks <- data.frame(rank_id = names(nmf_params$ranks),
#'   #                                   rank = nmf_params$ranks,
#'   #                                   Best_Frob_error = sapply(complete_eval, function(x) min(x$Frob_error)),
#'   #                                   Best_n_iter     = sapply(complete_eval, function(x) x$iterations[which.min(x$Frob_error)]),
#'   #                                   stringsAsFactors = FALSE)
#'
#'
#'   #frob_errors <- DataFrame(getFrobError(dec.matrix))
#'   #colnames(frob.errors) <- as.character(k.min:k.max)
#'
#'   #----------------------------------------------------------------------------#
#'   #                       Return nmfExperiment_lite object                     #
#'   #----------------------------------------------------------------------------#
#'   nmfExperiment_lite(input_matrix = input_matrix,
#'                      WMatrix      = lapply(complete_eval, "[[" , "W"),
#'                      HMatrix      = lapply(complete_eval, "[[" , "H"),
#'                      FrobError    = frob_errors,
#'                      OptKStats    = OptKStats,
#'                      OptK         = OptK)
#'   #FeatureStats              = "data.frame",
#'   #SignatureSpecificFeatures = "list" )
#'   #Factorization_ranks       = Factorization_ranks)
#' }
#'
#' environment(runNMFtensor_lite_eval) <- asNamespace("Bratwurst")
#' nmf_exp_eval <- runNMFtensor_lite_eval(norm_mat_list[[1]],
#'                              #ranks = 2:10,
#'                              ranks = c(2,4,6,8,10,12),
#'                              n_initializations     = 10,
#'                              iterations            = 10^4,
#'                              convergence_threshold = 5)
#' nmf_exp_eval
#' str(nmf_exp_eval$k6$W_eval)
#'
#'
#'
#' # nmf_exp_eval_OptKStats_NMF <- compute_OptKStats_NMF(nmf_exp_eval)
#' # nmf_exp_eval_OptKStats_NMF
#'
#'
#'
#'
#'
#' #'
#' #'
#' #'
#' #'
#' #'
#'
#'
#'
#'

#'
#'











