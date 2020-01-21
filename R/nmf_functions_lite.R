# Copyright © 2015-2020  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#==============================================================================#
#                         NMF GPU Wrapper - FUNCTIONS                          #
#==============================================================================#

#==============================================================================#
#                      NMF tensorflow Wrapper - FUNCTION                       #
#==============================================================================#

#' Single view NMF
#'
#' Computes NMF on tensorflow using the reticulate framework, uses a non-negative matrix as input
#'
#' @param X Input matrix, should be a non-negative matrix
#' @param n_initializations Number of initializations to evaluate
#' @param iterations Maximum number of iterations to run for every initialization
#' @param convergence_threshold The factorization stops, if the convergence test is constant for this number of iterations
#'
#' @return A nmfExperiment_lite object, containg the initial matrix X and the faactorized matrices W and H, along with factorization metrics
#'
#' @export
#'
#' @examples
#' nmf_exp <- runNMFtensor_lite(matrix(1:1000, ncol = 10),
#'                              ranks = 2:5,
#'                              n_initializations     = 10,
#'                              iterations            = 10^4,
#'                              convergence_threshold = 40)
#' nmf_exp
runNMFtensor_lite <- function (X,
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
  #print(frob_errors)
  #sapply(complete_eval, function(x) print(x$Frob_error))
  #print(sapply(complete_eval, function(x) x$Frob_error))
  # Factorization_ranks <- data.frame(rank_id = names(nmf_params$ranks),
  #                                   rank = nmf_params$ranks,
  #                                   Best_Frob_error = sapply(complete_eval, function(x) min(x$Frob_error)),
  #                                   Best_n_iter     = sapply(complete_eval, function(x) x$iterations[which.min(x$Frob_error)]),
  #                                   stringsAsFactors = FALSE)


  #frob_errors <- DataFrame(getFrobError(dec.matrix))
  #colnames(frob.errors) <- as.character(k.min:k.max)

  #----------------------------------------------------------------------------#
  #                       Return nmfExperiment_lite object                     #
  #----------------------------------------------------------------------------#
  nmfExperiment_lite(input_matrix = input_matrix,
                     WMatrix      = lapply(complete_eval, "[[" , "W"),
                     HMatrix      = lapply(complete_eval, "[[" , "H"),
                     FrobError    = frob_errors,
                     OptKStats    = "list",
                     OptK         = "numeric")
                     #FeatureStats              = "data.frame",
                     #SignatureSpecificFeatures = "list" )
                     #Factorization_ranks       = Factorization_ranks)
}

#==============================================================================#
#                               Getter FUNCTIONS                               #
#==============================================================================#
# getFrobError
# getHMatrixList
# getWMatrixList

#'
#' #==============================================================================#
#' #             Criteria for optimal factorization rank - FUNCTIONS              #
#' #==============================================================================#
#' #' Compute basic statistics for Frobenius Errors
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' computeFrobErrorStats <- function(nmf.exp) {
#'   frob.errorMatrix <- as.matrix(FrobError(nmf.exp))
#'   min.frobError <- apply(frob.errorMatrix, 2, function(x) min(x, na.rm = T))
#'   sd.frobError <- apply(frob.errorMatrix, 2, function(x) sd(x, na.rm = T))
#'   mean.frobError <- colMeans(frob.errorMatrix, na.rm = T)
#'   cv.frobError <- sd.frobError / mean.frobError
#'   frobError.data <- DataFrame("k" = as.numeric(names(min.frobError)),
#'                               "min" = min.frobError,
#'                               "mean" = mean.frobError,
#'                               "sd" = sd.frobError,
#'                               "cv" = cv.frobError)
#'   nmf.exp <- setOptKStats(nmf.exp, frobError.data)
#'   return(nmf.exp)
#' }
#'
#' #' Compute p-value with t-test for running K
#' #'
#' #' @param frobError.matrix
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' computePValue4FrobError <- function(frobError.matrix) {
#'   n.col <- ncol(frobError.matrix)
#'   p <- mapply(1:(n.col - 1), 2:n.col, FUN = function(i, j) {
#'     p <- t.test(frobError.matrix[, i], frobError.matrix[, j])
#'     return(p$p.value)
#'   })
#'   #p <- p.adjust(unlist(p), method = 'BH')
#'   return(-log10(unlist(p)))
#' }
#'
#' #' Cosine similarity
#' #'
#' #' @param a
#' #' @param b
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' cosineSim <- function(a, b){
#'   return(sum(a * b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ))
#' }
#'
#' #' Cosine distance
#' #'
#' #' @param a
#' #' @param b
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' cosineDist <- function(a, b){
#'   return(1 - cosineSim(a, b))
#' }
#'
#' #' Create distance matrix with cosine similarity
#' #'
#' #' @param in.matrix
#' #' @param in.dimension
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' cosineDiss <- function(in.matrix, in.dimension=2){
#'   if(in.dimension == 1) in.matrix <- t(in.matrix)
#'   cosineDist.list <- lapply(1:ncol(in.matrix), function(i.outCol) {
#'     cosine.dists <- unlist(lapply(1:ncol(in.matrix), function(i.innCol) {
#'       cosineDist(in.matrix[, i.outCol], in.matrix[, i.innCol])
#'     }))
#'   })
#'   diss.matrix <- do.call(rbind, cosineDist.list)
#'   return(round(diss.matrix, digits = 14))
#' }
#'
#' #' Create distance matrix with cosine similarity with matrix operations
#' #'
#' #' @param in.matrix
#' #' @param in.dimension
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' cosineDissMat <- function(in.matrix, in.dimension=2){
#'   if(in.dimension == 1) in.matrix <- t(in.matrix)
#'   squaredVectorSum <- apply(in.matrix, 2, function(m) { sqrt(sum(m * m)) })
#'   squaredVectorProduct <- squaredVectorSum %*% t(squaredVectorSum)
#'   squaredInputSum <- t(in.matrix) %*% in.matrix
#'   # sum(a*b) for any a,b in M
#'   diss.matrix <- 1 - squaredInputSum / squaredVectorProduct
#'   # CosineDistance = 1 - CosineSimilarity
#'   return(round(diss.matrix, digits = 14))
#' }
#'
#' #' Compute Alexandrov Criterion --> Silhoutte Width
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return
#' #'
#' #' @importFrom cluster pam
#' #' @export
#' #'
#' #' @examples
#' computeSilhoutteWidth <- function(nmf.exp) {
#'   sil.vec <- lapply(WMatrixList(nmf.exp), function(WMatrix.list) {
#'     nan.bool <- sapply(WMatrix.list, function(m) !any(is.nan(m)))
#'     concat.matrix <- do.call(cbind, WMatrix.list[which(nan.bool)])
#'     dist.matrix <- cosineDissMat(as.matrix(concat.matrix))
#'     my.pam <- pam(dist.matrix, k = ncol(WMatrix.list[[1]]),  diss = T)
#'     sil.sum <- sum(my.pam$silinfo$widths[, "sil_width"])
#'     sil.mean <- mean(my.pam$silinfo$widths[, "sil_width"])
#'     return(DataFrame(sumSilWidth = sil.sum,
#'                      meanSilWidth = sil.mean))
#'   })
#'   sil.vec <- do.call(rbind, sil.vec)
#'   nmf.exp <- setOptKStats(nmf.exp, cbind(OptKStats(nmf.exp), sil.vec))
#'   return(nmf.exp)
#' }
#'
#' #' Compute Cophenetic correlation coefficient, TO BE IMPROVED
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' computeCopheneticCoeff <- function(nmf.exp) {
#'   coph.coeff <- lapply(WMatrixList(nmf.exp), function(WMatrix.list) {
#'     nan.bool <- sapply(WMatrix.list, function(m) !any(is.nan(m)))
#'     concat.matrix <- do.call(cbind, WMatrix.list[which(nan.bool)])
#'     #concat.matrix <- do.call(cbind, WMatrix.list)
#'     dist.matrix <- cosineDissMat(as.matrix(concat.matrix))
#'     my.hclust <- hclust(as.dist(dist.matrix))
#'     dist.cophenetic <- as.matrix(cophenetic(my.hclust))
#'     # take distance matrices without diagonal elements
#'     diag(dist.matrix) <- NA
#'     dist.matrix <- dist.matrix[which(!is.na(dist.matrix))]
#'     diag(dist.cophenetic) <- NA
#'     dist.cophenetic <- dist.cophenetic[which(!is.na(dist.cophenetic))]
#'     return(cor(cbind(dist.cophenetic, dist.matrix))[1, 2])
#'   })
#'   coph.coeff <- DataFrame(copheneticCoeff = unlist(coph.coeff))
#'   nmf.exp <- setOptKStats(nmf.exp, cbind(OptKStats(nmf.exp), coph.coeff))
#'   return(nmf.exp)
#' }
#'
#' #' Compute amari type distance between two matrices
#' #'
#' #' @param matrix.A,matrix.B of the same dimensionality
#' #'
#' #' @return The amari type distance of matrix.A & matrix.B according
#' #'        to [Wu et. al, PNAS 2016]
#' #'
#' #' @references \url{http://www.pnas.org/content/113/16/4290.long}
#' #'
#' #' @export
#' #'
#' #' @examples
#' amariDistance <- function(matrix.A, matrix.B) {
#'   K <- dim(matrix.A)[2]
#'   C <- cor(matrix.A, matrix.B)
#'   return(1 - (sum(apply(C, 1, max)) + sum(apply(C, 2, max))) / (2 * K))
#' }
#'
#' #' Compute Amari Distances from [Wu et. al, PNAS 2016]
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return The average Amari-type error for each k
#' #'
#' #' @references \url{http://www.pnas.org/content/113/16/4290.long}
#' #'
#' #' @export
#' #'
#' #' @examples
#' computeAmariDistances <- function(nmf.exp){
#'   distance.averages <- lapply(WMatrixList(nmf.exp), function(matrices) {
#'     B <- length(matrices)
#'     distances.list <- unlist(lapply(1:(B - 1), function(b) {
#'       distances <- lapply((b + 1):B, function(b.hat) {
#'         amariDistance(matrices[[b]], matrices[[b.hat]])
#'       })
#'     }))
#'     return(mean(distances.list[!is.na(distances.list)]))
#'     # is.na to exclude corrupted matrices
#'   })
#'   distance.averages <- DataFrame(meanAmariDist = unlist(distance.averages))
#'   nmf.exp <- setOptKStats(nmf.exp, cbind(OptKStats(nmf.exp), distance.averages))
#'   return(nmf.exp)
#' }
#'
#' local.minima <- function(x)
#'   ifelse(dplyr::lag(x) >= x & dplyr::lead(x) >= x, TRUE, FALSE)
#' local.maxima <- function(x)
#'   ifelse(dplyr::lag(x) <= x & dplyr::lead(x) <= x, TRUE, FALSE)
#'
#'
#' #' Title
#' #'
#' #' @param nmf.exp
#' #' @param verbose
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' proposeOptK <- function(nmf.exp, verbose = FALSE){
#'   OptKStats_DF <- OptKStats(nmf.exp)
#'   if(ncol(OptKStats_DF) & nrow(OptKStats_DF)){
#'     indCopheneticCoeff <- which(local.maxima(OptKStats_DF$copheneticCoeff))
#'     indMeanAmariDist <- which(local.minima(OptKStats_DF$meanAmariDist))
#'     return(list(
#'       intersectK = OptKStats_DF$k[intersect(indCopheneticCoeff,
#'                                             indMeanAmariDist)],
#'       unionK = OptKStats_DF$k[union(indCopheneticCoeff, indMeanAmariDist)]))
#'   } else {
#'     if(verbose) cat("proposeOptK::warning:OptKStats not yet computed.\n")
#'     return(NULL)
#'   }
#' }
#'
#' #==============================================================================#
#' #                         H-MATRIX ANALYSIS FUNCTIONS                          #
#' #==============================================================================#
#' # Compute unbaised signature names by applying a row k-means
#' # and classify signature according to cluster and highest cluster mean.
#' #' Compute unsupervised signature names from colData and H-Matrix for given K
#' #'
#' #' @param nmf.exp
#' #' @param k.opt
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' getSignatureNames <- function(nmf.exp, k.opt) {
#'   H <- HMatrix(nmf.exp, k = k.opt)
#'   sig.names <- lapply(1:nrow(H), function(i) {
#'     k.cluster <- kmeans(as.numeric(H[i, ]), 2)
#'     n.kCluster <- which(k.cluster$centers == max(k.cluster$centers))
#'     samples.kCluster <- which(k.cluster$cluster == n.kCluster)
#'     meta.data <- colData(nmf.exp)[samples.kCluster, ]
#'     n.col <- ncol(meta.data)
#'     if (n.col > 1) {
#'       sigName.combs <- apply(as.data.frame(meta.data), 1, function(x){
#'         paste(x[2:n.col], collapse = " ")
#'       })
#'       sigName.combs <- sort(table(sigName.combs), decreasing = T)
#'       sig.name <- names(sigName.combs)[1]
#'       # TODO: Might be improvable by using exposures from H-Matrix.
#'       # Exposure proportion computation!?
#'       if(length(sigName.combs) > 1) {
#'         sig.prop <- round(sigName.combs / sum(sigName.combs), 2)
#'         sig.prop <- paste(names(sig.prop), sig.prop)
#'         sig.name <- paste(sig.prop, collapse = "\n")
#'       }
#'     } else {
#'       sig.name <- sprintf("Signature %s", i)
#'     }
#'     return(sig.name)
#'   })
#'   return(unlist(sig.names))
#' }
#'
#'
#' #==============================================================================#
#' #                         W-MATRIX ANALYSIS FUNCTIONS:                         #
#' #                             FEATURE SELECTION                                #
#' #==============================================================================#
#' #' Compute 'shannon' entropy per region.
#' #' High Entropy means highly specific for one signature.
#' #'
#' #' @param matrix
#' #'
#' #' @return
#' #'
#' #' @examples
#' computeEntropy <- function(matrix) {
#'   matrix.relativ <- t(apply(matrix, 1, function(x) x / sum(x)))
#'   matrix.entropy <- apply(matrix.relativ, 1, function(x) {
#'     p <- x * log2(length(x) * x)
#'     p[is.na(p)] <- 0
#'     h <- sum(p)
#'     return(h)
#'   })
#'   return(matrix.entropy)
#' }
#'
#' #' Computes entropy for each feature in optimal K W-matrix
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' computeEntropy4OptK <- function(nmf.exp) {
#'   W <- WMatrix(nmf.exp, OptK(nmf.exp))
#'   W.entropy <- computeEntropy(W)
#'   nmf.exp <- setFeatureStats(nmf.exp, DataFrame("entropy" = W.entropy))
#'   return(nmf.exp)
#' }
#'
#' #' Compute delta between each column (signature) per row
#' #'
#' #' @param matrix
#' #'
#' #' @return
#' #'
#' #' @examples
#' computeAbsDelta <- function(matrix) {
#'   delta.regions <- lapply(1:ncol(matrix), function(k) {
#'     delta.vec <- matrix[, k] - (rowSums(matrix[, -k]))
#'     return(delta.vec)
#'   })
#'   delta.regions <- do.call(cbind, delta.regions)
#'   return(delta.regions)
#' }
#'
#' #' Computes absolut delta per feature for each signature given optimal K
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' computeAbsDelta4OptK <- function(nmf.exp) {
#'   W <- WMatrix(nmf.exp, OptK(nmf.exp))
#'   W.absDelta <- computeAbsDelta(W)
#'   nmf.exp <- setFeatureStats(nmf.exp, DataFrame("absDelta" = W.absDelta))
#'   return(nmf.exp)
#' }
#'
#' #' Compute the coeficient of variation per row in a matrix.
#' #'
#' #' @param matrix
#' #'
#' #' @return
#' #'
#' #' @examples
#' computeCoefVar <- function(matrix) {
#'   apply(matrix, 1, function(r) sd(r) / mean(r))
#' }
#'
#' #' Perform Kmeans on rows of a matrix to classify them into column (singature) combinations
#' #'
#' #' @param matrix
#' #'
#' #' @return
#' #'
#' #' @importFrom cluster silhouette
#' #' @export
#' #'
#' #' @examples
#' performRowKmeans <- function(matrix) {
#'   k.row <- apply(matrix, 1, function(x) {
#'     # Perform sigmoidal transformation to achieve better clustering
#'     x.trans <- sigmoidTransform(x)
#'     k.cluster <- kmeans(x.trans, 2)
#'     d <- dist(as.data.frame(x.trans))
#'     sil.mean <- mean(silhouette(k.cluster$cluster, d)[, 3])
#'     cluster.deltaMean <- mean(x[k.cluster$cluster == 1]) -
#'       mean(x[k.cluster$cluster == 2])
#'     return(list(centers = t(k.cluster$centers),
#'                 silmean = sil.mean,
#'                 explainedVar = k.cluster$betweenss / k.cluster$totss,
#'                 oddsVar = sum(k.cluster$withinss) / k.cluster$betweenss,
#'                 attribution = k.cluster$cluster,
#'                 deltaMean = cluster.deltaMean))
#'   })
#'   return(k.row)
#' }
#'
#' #' Compute Signature Combinations for W-matrix given optimal K
#' #'
#' #' @param nmf.exp
#' #' @param var.thres
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' computeFeatureStats <- function(nmf.exp, var.thres = 0.25) {
#'   W <- WMatrix(nmf.exp, OptK(nmf.exp))
#'   # Determine Region contributions
#'   k.row <- performRowKmeans(W)
#'   # Get kmeans stats
#'   k.explainedVar <- unlist(lapply(k.row, function(r) r$explainedVar))
#'   k.oddsVar <- unlist(lapply(k.row, function(r) r$oddsVar))
#'   # Extract mean silhouette width for each row.
#'   k.silmean <- unlist(lapply(k.row, function(r) r$silmean))
#'   # Compute delta between cluster centers (transformed data!)
#'   k.deltaCenter <-
#'     unlist(lapply(k.row, function(r) r$centers[1] - r$centers[2]))
#'   # Extract difference of cluster means
#'   k.deltaMean <- unlist(lapply(k.row, function(r) r$deltaMean))
#'   # Extract Signature combinations generated by k-means
#'   k.attribution <- lapply(k.row, function(r) abs(r$attribution))
#'   k.attribution <- do.call(rbind, k.attribution)
#'   k.ids <- apply(k.attribution, 1, function(r) paste(r, collapse = ""))
#'   # Determine regions which contribute to all signatures.
#'   k.varCoef <- computeCoefVar(W)
#'   all.signature <- which(k.varCoef < var.thres)
#'   k.ids[all.signature] <- gsub("2", "1", k.ids[all.signature])
#'   # Set FeatureStats
#'   feature.stats <- DataFrame("cluster" = k.ids,
#'                              "deltaCenters" = k.deltaCenter,
#'                              "deltaMean" = k.deltaMean,
#'                              "explainedVar" = k.explainedVar,
#'                              "oddsVar" = k.oddsVar,
#'                              "coefVar" = k.varCoef,
#'                              "meanSil" = k.silmean)
#'   # Map reverse cluster definitions to each other,
#'   # including sign change for delta means
#'   ids <- sort(unique(k.ids))[-1] ## this is not robust!!
#'   ## must be replaced by
#'   # equalPattern <- paste(rep("1", ncol(W)), collapse = "")
#'   # ids <- setdiff(sort(unique(k.ids)), equalPattern)
#'   ids1 <- ids[1:(length(ids) / 2)]
#'   ids2 <-  gsub("0", "2", gsub("2", "1", gsub("1", "0", ids1)))
#'   conv.id <- data.frame("id1" = ids1, "id2" = ids2)
#'   i.conv <- match(feature.stats$cluster, conv.id$id2)
#'   feature.stats$cluster[!is.na(i.conv)] <-
#'     as.character(conv.id$id1[i.conv[!is.na(i.conv)]])
#'   feature.stats$deltaMean[!is.na(i.conv)] <-
#'     (-1) * feature.stats$deltaMean[!is.na(i.conv)]
#'   feature.stats$deltaCenters[!is.na(i.conv)] <-
#'     (-1) * feature.stats$deltaCenters[!is.na(i.conv)]
#'   # Re-write cluster ids in a more useful binary code
#'   i <- which(feature.stats$deltaCenters > 0 &
#'                grepl(feature.stats$cluster, pattern = "2"))
#'   feature.stats$cluster[i] <- gsub("2", "0", feature.stats$cluster[i])
#'   i <- which(feature.stats$deltaCenters < 0 &
#'                grepl(feature.stats$cluster, pattern = "2"))
#'   feature.stats$cluster[i] <-
#'     gsub("2", "1", gsub("1", "0", feature.stats$cluster[i]))
#'   # Add featureStats to NMF experiment
#'   nmf.exp <- setFeatureStats(nmf.exp, feature.stats)
#'   return(nmf.exp)
#' }
#'
#' #' Title
#' #'
#' #' @param k.ids
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' getIndex4SignatureCombs <- function(k.ids) {
#'   # Get for all signature combinations corresponding regions.
#'   ids <- sort(unique(k.ids))[-1]
#'   n.ids <- length(ids) + 1
#'   i.regions <- lapply(1:(n.ids / 2), function(i) {
#'     sub.id <- ids[c(i, n.ids - i)]
#'     i.region <- which(k.ids %in% sub.id)
#'     return(i.region)
#'   })
#'   names(i.regions) <- ids[1:(n.ids / 2)]
#'   # Get regions for all signature id.
#'   allSig.id <- sort(unique(k.ids))[1]
#'   i.regions[[allSig.id]] <- which(k.ids %in% allSig.id)
#'   return(i.regions)
#' }
#'
#' #' Title
#' #'
#' #' @param W
#' #' @param i.regions
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' computeMeanDiff4SignatureCombs <- function(W, i.regions) {
#'   # Compute Row mean diff between Cluster 1 and 2.
#'   w.diffs <- lapply(1:length(i.regions), function(i) {
#'     w <- W[i.regions[[i]], ]
#'     if (!grepl(names(i.regions)[i], pattern = "2")) {
#'       w.diff <- rowMeans(w)
#'     } else {
#'       signature.comb <-
#'         as.numeric(unlist(strsplit(names(i.regions)[[i]], split = "")))
#'       w.mean1 <- rowMeans(w[, which(signature.comb == 1)])
#'       w.mean2 <- rowMeans(w[, which(signature.comb == 2)])
#'       w.diff <- w.mean1 - w.mean2
#'     }
#'     return(w.diff)
#'   })
#'   names(w.diffs) <- names(i.regions)
#' }
#'
#' ### Get number of peaks per signature combination and create a barplot.
#' #' Title
#' #'
#' #' @param i.regions
#' #' @param w.diffs
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' getSignatureCombCounts <- function(i.regions, w.diffs) {
#'   n.peaks <- lapply(names(i.regions), function(region.n) {
#'     n.peak1 <- sum(w.diffs[[region.n]] > 0)
#'     n.peak2 <- sum(w.diffs[[region.n]] < 0)
#'     n.peak <- c(n.peak1, n.peak2)
#'     n.peak <- melt(n.peak)
#'     n.peak$clusterId <- c("1", "2")
#'     n.peak$sigCombId <- region.n
#'     return(n.peak)
#'   })
#'   n.peaks <- do.call(rbind, n.peaks)
#'   return(n.peaks)
#' }
#'
#'
#' #' Normalize the signatures matrix (W)
#' #'
#' #' After column normalization of the matrix W, the inverse factors are
#' #' mutiplied with the rows of H in order to keep the matrix product W*H
#' #' constant.
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return A data structure of type nmfExperiment
#' #'
#' #' @importFrom YAPSA normalize_df_per_dim
#' #' @export
#' #'
#' #' @examples
#' #'  NULL
#' #'
#' normalizeW <- function(nmf.exp){
#'   # account for WMatrixList and HMatrixList
#'   all_list <- lapply(seq_along(WMatrixList(nmf.exp)), function(k_ind){
#'     k_list <-
#'       lapply(seq_along(WMatrixList(nmf.exp)[[k_ind]]), function(init_ind){
#'         tempW <- WMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         tempH <- HMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         normFactor <- colSums(tempW)
#'         # catch errors associated with NaNs in W or H
#'         if (any(is.nan(normFactor))){
#'           return(list(W = tempW,
#'                       H = tempH))
#'         }else{
#'           newSigs <- as.matrix(normalize_df_per_dim(tempW, 2))
#'           newExpo <- tempH * normFactor
#'           #newV <- newSigs %*% newExpo
#'           #oldV <- tempW %*% tempH
#'           return(list(W = newSigs,
#'                       H = newExpo))
#'         }
#'       })
#'     names(k_list) <- names(WMatrixList(nmf.exp)[[k_ind]])
#'     return(k_list)
#'   })
#'   names(all_list) <- names(WMatrixList(nmf.exp))
#'   thisWMatrixList <- lapply(all_list, function(current_k_list){
#'     kWMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$W)
#'     })
#'   })
#'   nmf.exp <- setWMatrixList(nmf.exp, thisWMatrixList)
#'   thisHMatrixList <- lapply(all_list, function(current_k_list){
#'     kHMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$H)
#'     })
#'   })
#'   nmf.exp <- setHMatrixList(nmf.exp, thisHMatrixList)
#'   return(nmf.exp)
#' }
#'
#' #' Normalize the signatures matrix (H)
#' #'
#' #' After row normalization of the matrix H, the inverse factors are
#' #' mutiplied with the columns of W in order to keep the matrix product W*H
#' #' constant.
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return A data structure of type nmfExperiment
#' #'
#' #' @importFrom YAPSA normalize_df_per_dim
#' #' @export
#' #'
#' #' @examples
#' #'  NULL
#' #'
#' normalizeH <- function(nmf.exp){
#'   # account for WMatrixList and HMatrixList
#'   all_list <- lapply(seq_along(WMatrixList(nmf.exp)), function(k_ind){
#'     k_list <-
#'       lapply(seq_along(WMatrixList(nmf.exp)[[k_ind]]), function(init_ind){
#'         tempW <- WMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         tempH <- HMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         normFactor <- rowSums(tempH)
#'         newExpo <- as.matrix(normalize_df_per_dim(tempH, 1))
#'         newSigs <- tempW * normFactor
#'         return(list(W = newSigs,
#'                     H = newExpo))
#'       })
#'     names(k_list) <- names(WMatrixList(nmf.exp)[[k_ind]])
#'     return(k_list)
#'   })
#'   names(all_list) <- names(WMatrixList(nmf.exp))
#'   thisWMatrixList <- lapply(all_list, function(current_k_list){
#'     kWMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$W)
#'     })
#'   })
#'   nmf.exp <- setWMatrixList(nmf.exp, thisWMatrixList)
#'   thisHMatrixList <- lapply(all_list, function(current_k_list){
#'     kHMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$H)
#'     })
#'   })
#'   nmf.exp <- setHMatrixList(nmf.exp, thisHMatrixList)
#'   return(nmf.exp)
#' }
#'
#' #' Regularize the signatures matrix (H)
#' #'
#' #' After row regularization of the matrix H, the inverse factors are
#' #' mutiplied with the columns of W in order to keep the matrix product W*H
#' #' constant.
#' #'
#' #' @param nmf.exp
#' #'
#' #' @return A data structure of type nmfExperiment
#' #'
#' #' @export
#' #'
#' #' @examples
#' #'  NULL
#' #'
#' regularizeH <- function(nmf.exp){
#'   # account for WMatrixList and HMatrixList
#'   all_list <- lapply(seq_along(WMatrixList(nmf.exp)), function(k_ind){
#'     k_list <-
#'       lapply(seq_along(WMatrixList(nmf.exp)[[k_ind]]), function(init_ind){
#'         tempW <- WMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         tempH <- HMatrixList(nmf.exp)[[k_ind]][[init_ind]]
#'         normFactor <- rowMax(tempH)
#'         newExpo <- tempH / normFactor
#'         newSigs <- tempW * normFactor
#'         return(list(W = newSigs,
#'                     H = newExpo))
#'       })
#'     names(k_list) <- names(WMatrixList(nmf.exp)[[k_ind]])
#'     return(k_list)
#'   })
#'   names(all_list) <- names(WMatrixList(nmf.exp))
#'   thisWMatrixList <- lapply(all_list, function(current_k_list){
#'     kWMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$W)
#'     })
#'   })
#'   nmf.exp <- setWMatrixList(nmf.exp, thisWMatrixList)
#'   thisHMatrixList <- lapply(all_list, function(current_k_list){
#'     kHMatrixList <- lapply(current_k_list, function(current_entry){
#'       return(current_entry$H)
#'     })
#'   })
#'   nmf.exp <- setHMatrixList(nmf.exp, thisHMatrixList)
#'   return(nmf.exp)
#' }
#'
#'
#' #' Compute signature specific features
#' #'
#' #' @param nmf.exp A NMF experiment object
#' #' @param rowDataId The index of the rowData(nmf.exp) data.frame that should be used for
#' #'  feature extraction. In case rowData(nmf.exp)[,1] is a GRanges or a related object like
#' #'  GenomicInteractions this parameter can be ignored
#' #'
#' #' @return nmf.exp with filles SignatureSpecificFeatures container
#' #'
#' #'
#' #' @export
#' #'
#' #' @examples
#' #'
#' computeSignatureSpecificFeatures <- function(nmf.exp, rowDataId = 3){
#'   if (length(OptK(nmf.exp)) == 0){
#'     stop("You need to first define an optimal k before being able to compute
#'            signature specific features!")
#'   } else {
#'     if (nrow(FeatureStats(nmf.exp)) == 0){
#'       message("Computing feature stats...")
#'       nmf.exp <- computeFeatureStats(nmf.exp)
#'     }
#'     fstats <- FeatureStats(nmf.exp)
#'     # identify unique cluster membership strings
#'     clusterMemberships <-
#'       sapply(unique(as.character(fstats$cluster)),
#'              function (x) lengths(regmatches(x, gregexpr("1", x))))
#'     sigSpecClusters <-
#'       sort(names(clusterMemberships[which(clusterMemberships == 1)]),
#'            decreasing = TRUE)
#'
#'     if (class(rowData(nmf.exp)[, 1]) %in% c("GRanges", "GInteractions",
#'                                             "GenomicInteractions")){
#'       signatureSpecificFeatures <- lapply(sigSpecClusters, function(ssc){
#'         features <- rowData(nmf.exp)[, 1][which(fstats$cluster == ssc)]
#'         return(features)
#'       })
#'     }else{
#'       signatureSpecificFeatures <- lapply(sigSpecClusters, function(ssc){
#'         features <- rowData(nmf.exp)[, rowDataId][which(fstats$cluster == ssc)]
#'         return(features)
#'       })
#'     }
#'     names(signatureSpecificFeatures) <- sigSpecClusters
#'     nmf.exp@SignatureSpecificFeatures <- signatureSpecificFeatures
#'     return(nmf.exp)
#'   }
#' }
#'
#'
