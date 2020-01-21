runNMFtensor_lite_eval <- function (X,
                               ranks,
                               n_initializations     = 10,
                               iterations            = 10^4,
                               convergence_threshold = 40){

  #----------------------------------------------------------------------------#
  #                                    Setup                                   #
  #----------------------------------------------------------------------------#
  # Check  data
  if (min(X) < 0 ) {
    stop("\nNegative values present in input matrix\n
         only non-negative matrices supported\n")
  }
  # Convert params to integer
  nmf_params <- lapply(list(ranks                 = ranks,
                            n_initializations     = n_initializations,
                            iterations            = iterations,
                            convergence_threshold = convergence_threshold),
                       as.integer)

  #----------------------------------------------------------------------------#
  #                Run single view NMF on tensorflow                           #
  #----------------------------------------------------------------------------#
  # Source NMF tensorflow python script
  source_python(file.path(system.file(package = "Bratwurst"), "python/nmf_tensor_lite.py"))

  # Run NMF
  #X <- input_matrix
  names(nmf_params$ranks) <- paste0("k", nmf_params$ranks)
  complete_eval <- lapply(nmf_params$ranks, function(k) {
    #k <- as.integer(k)

    print(Sys.time())
    cat("Factorization rank: ", k, "\n")

    k_eval <- NMF_tensor_py(matrix            = X,
                            rank              = k,
                            n_initializations = nmf_params$n_initializations,
                            iterations        = nmf_params$iterations,
                            stop_threshold    = nmf_params$convergence_threshold)
    names(k_eval) <- c("W", "H", "iterations", "Frob_error", "W_eval")
    k_eval$iterations <- unlist(k_eval$iterations)
    k_eval$Frob_error <- unlist(k_eval$Frob_error)


    print(paste("NMF converged after ", paste(k_eval$iterations, collapse = ","), "iterations"))
    return(k_eval)
  })
  return(complete_eval)

  #----------------------------------------------------------------------------#
  #                        Build NMF object slots                              #
  #----------------------------------------------------------------------------#
  # input matrix info
  input_matrix <- list(hash     = digest::digest(X),
                       dim      = dim(X),
                       colnames = colnames(X),
                       rownames = rownames(X))

  #names(dec.matrix) <- k.min:k.max
  frob_errors <- as.data.frame(do.call(cbind, lapply(complete_eval, "[[" , "Frob_error")))
  #return(frob_errors)
  #print(frob_errors)
  #sapply(complete_eval, function(x) print(x$Frob_error))
  #print(sapply(complete_eval, function(x) x$Frob_error))
  Factorization_ranks <- data.frame(rank_id = names(nmf_params$ranks),
                                    rank = nmf_params$ranks,
                                    Best_Frob_error = sapply(complete_eval, function(x) min(x$Frob_error)),
                                    Best_n_iter     = sapply(complete_eval, function(x) x$iterations[which.min(x$Frob_error)]),
                                    stringsAsFactors = FALSE)


  #frob_errors <- DataFrame(getFrobError(dec.matrix))
  #colnames(frob.errors) <- as.character(k.min:k.max)

  #----------------------------------------------------------------------------#
  #                       Return nmfExperiment_lite object                     #
  #----------------------------------------------------------------------------#
  nmfExperiment_lite(input_matrix = input_matrix,
                     WMatrix      = lapply(complete_eval, "[[" , "W"),
                     HMatrix      = lapply(complete_eval, "[[" , "H"),
                     FrobError    = frob_errors,
                     #OptKStats    = "list",
                     #OptK         = "numeric",
                     #FeatureStats              = "data.frame",
                     #SignatureSpecificFeatures = "list" )
                     Factorization_ranks       = Factorization_ranks)

}
environment(runNMFtensor_lite_eval) <- asNamespace("Bratwurst")
nmf_exp_eval <- runNMFtensor_lite_eval(norm_mat_list[[1]],
                             #ranks = 2,
                             ranks = c(2,4,6),
                             n_initializations     = 2,
                             iterations            = 10^4,
                             convergence_threshold = 5)
nmf_exp_eval
str(nmf_exp_eval$k6$W_eval)


compute_OptKStats_NMF <- function(complete_eval) {
  #----------------------------------------------------------------------------#
  #                            Frobenius error stats                           #
  #----------------------------------------------------------------------------#
  frob_errors_df <- as.data.frame(do.call(cbind, lapply(complete_eval, "[[" , "Frob_error")))
  min_frobError  <- apply(frob_errors_df, 2, function(x) min(x, na.rm = TRUE))
  sd_frobError   <- apply(frob_errors_df, 2, function(x) sd(x, na.rm = TRUE))
  mean_frobError <- colMeans(frob_errors_df, na.rm = TRUE)
  cv_frobError   <- sd_frobError / mean_frobError

  #----------------------------------------------------------------------------#
  #     compute Silhouette Width, Cophenetic Coeff and Amari Distances         #
  #----------------------------------------------------------------------------#
  sil_vec <- lapply(lapply(complete_eval, "[[", "W_eval"), function(WMatrix_list) {
    B <- length(WMatrix_list)

    # nan.bool <- sapply(WMatrix_list, function(m) !any(is.nan(m)))
    # concat_matrix <- do.call(cbind, WMatrix_list[which(nan.bool)])
    concat_matrix <- do.call(cbind, WMatrix_list)
    dist_matrix   <- cosineDissMat(as.matrix(concat_matrix))

    # compute Silhouette Width
    if (length(WMatrix_list) > 1) {
      #------------------------------------------------------------------------#
      #                         compute Silhouette Width                       #
      #------------------------------------------------------------------------#
      my_pam   <- cluster::pam(dist_matrix, k = ncol(WMatrix_list[[1]]),  diss = TRUE)
      sumSilWidth  <- sum(my_pam$silinfo$widths[, "sil_width"])
      meanSilWidth <- mean(my_pam$silinfo$widths[, "sil_width"])
      #------------------------------------------------------------------------#
      #                         compute Cophenetic Coeff                       #
      #------------------------------------------------------------------------#
      # compute Cophenetic Coeff
      my_hclust <- hclust(as.dist(dist_matrix))
      dist_cophenetic <- as.matrix(cophenetic(my_hclust))
      # take distance matrices without diagonal elements
      diag(dist_matrix) <- NA
      dist_matrix <- dist_matrix[which(!is.na(dist_matrix))]
      diag(dist_cophenetic) <- NA
      dist_cophenetic <- dist_cophenetic[which(!is.na(dist_cophenetic))]
      copheneticCoeff = unlist(cor(cbind(dist_cophenetic, dist_matrix))[1, 2])

      #------------------------------------------------------------------------#
      #                         compute Amari Distances                        #
      #------------------------------------------------------------------------#
      distances_list <- unlist(lapply(1:(B - 1), function(b) {
        distances <- lapply((b + 1):B, function(b.hat) {
          amariDistance(WMatrix_list[[b]], WMatrix_list[[b.hat]])
        })
      }))
      meanAmariDist <- unlist(mean(distances_list))

    } else {
      sumSilWidth     <- NA
      meanSilWidth    <- NA
      copheneticCoeff <- NA
      meanAmariDist   <- NA
    }
    return(data.frame(sumSilWidth     = sumSilWidth,
                      meanSilWidth    = meanSilWidth,
                      copheneticCoeff = copheneticCoeff,
                      meanAmariDist   = meanAmariDist))
  })
  sil_vec <- dplyr::bind_rows(sil_vec)

  #----------------------------------------------------------------------------#
  #                         compute Amari Distances                            #
  #----------------------------------------------------------------------------#

  #----------------------------------------------------------------------------#
  #                         Return optimal K stats                             #
  #----------------------------------------------------------------------------#
  data.frame(k   = as.numeric(sub("^k", "", colnames(frob_errors_df))),
            min  = min_frobError,
            mean = mean_frobError,
            sd   = sd_frobError,
            cv   = cv_frobError,
            sumSilWidth     = sil_vec$sumSilWidth,
            meanSilWidth    = sil_vec$meanSilWidth,
            copheneticCoeff = sil_vec$copheneticCoeff,
            meanAmariDist   = sil_vec$meanAmariDist)
}
nmf_exp_eval_OptKStats_NMF <- compute_OptKStats_NMF(nmf_exp_eval)
nmf_exp_eval_OptKStats_NMF

library(tidyr)




local.minima <- function(x)
  ifelse(dplyr::lag(x) >= x & dplyr::lead(x) >= x, TRUE, FALSE)
local.maxima <- function(x)
  ifelse(dplyr::lag(x) <= x & dplyr::lead(x) <= x, TRUE, FALSE)


#' Title
#'
#' @param nmf.exp
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
proposeOptK <- function(nmf.exp, verbose = FALSE){
  OptKStats_DF <- OptKStats(nmf.exp)
  if(ncol(OptKStats_DF) & nrow(OptKStats_DF)){
    indCopheneticCoeff <- which(local.maxima(OptKStats_DF$copheneticCoeff))
    indMeanAmariDist <- which(local.minima(OptKStats_DF$meanAmariDist))
    return(list(
      intersectK = OptKStats_DF$k[intersect(indCopheneticCoeff,
                                            indMeanAmariDist)],
      unionK = OptKStats_DF$k[union(indCopheneticCoeff, indMeanAmariDist)]))
  } else {
    if(verbose) cat("proposeOptK::warning:OptKStats not yet computed.\n")
    return(NULL)
  }
}





gg_plotKStats <- function(nmf.exp,
                       plot.vars = c("FrobError", "cv", "sumSilWidth",
                                     "meanSilWidth", "copheneticCoeff",
                                     "meanAmariDist"),
                       optK = NULL) {
  frobError.df <- melt(as.data.frame(FrobError(nmf.exp)))
  frobError.df[, 1] <- factor(gsub("X", "", frobError.df[, 1]))
  frobError.df <- data.frame("k" = frobError.df[, 1],
                             "variable" = "FrobError",
                             "value" = frobError.df[, 2])
  optKStats.df <- melt(as.data.frame(OptKStats(nmf.exp)), id.vars = "k")
  meanError.df <- optKStats.df[optKStats.df$variable == "mean", ]
  meanError.df$variable <- unique(frobError.df$variable)
  #plot.vars <- as.character(unique(optKStats.df$variable)[-1:-3])
  optKStats.df <- optKStats.df[optKStats.df$variable %in% plot.vars, ]
  optKStats.df <- rbind(frobError.df, optKStats.df)
  optKStats.df$k <- as.numeric(as.character(optKStats.df$k))
  gg.optK <- ggplot() + geom_point(data = optKStats.df, aes(x = k, y = value),
                                   col = "black", size = 0.75)
  gg.optK <- gg.optK + geom_point(data = meanError.df,
                                  aes(x = k, y = value), col = "red")
  gg.optK <- gg.optK + facet_wrap(~variable, scales = "free_y")
  gg.optK <- gg.optK + xlab("K") + ylab("") + theme_bw() + science_theme()
  gg.optK <- gg.optK + theme(strip.background = element_rect(fill = "white"))
  if(!is.null(optK)){
    if(min(optKStats.df$k) <= optK & optK <= max(optKStats.df$k)){
      rect_df <- data.frame(xmin = optK - 0.5, xmax = optK + 0.5,
                            ymin = -Inf, ymax = Inf)
      gg.optK <- gg.optK +
        geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                      ymax = ymax),
                  fill = "grey50", alpha = 0.3, inherit.aes = FALSE)
    }
  }
  return(gg.optK)
}








