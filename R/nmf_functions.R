# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#==============================================================================#
#                         NMF GPU Wrapper - FUNCTIONS                          #
#==============================================================================#
#' Title
#'
#' @param matrix
#' @param tmp.path
#' @param sep
#'
#' @return
#'
#' @importFrom RcppCNPy npySave
#'
#' @examples
writeTmpMatrix <- function(matrix, tmp.path = "/tmp/nmf_tmp",
                           sep = " ", binary=FALSE) {
  # Write matrix to tmp file.
  matrix <- round(matrix, 7)
  dir.create(tmp.path)
  # Rmv final slash from path --> will lead to an error of nmfGpu wrapper
  tmp.path <- sub("/$", "", tmp.path)

  if(!binary){
    tmpMatrix.path <- tempfile(tmpdir = tmp.path, fileext = ".txt")
    write.table(matrix, file = tmpMatrix.path, quote = F,
                col.names = F, row.names = F, sep = sep)
  } else {
    tmpMatrix.path <- tempfile(tmpdir = tmp.path, fileext = ".npy")
    npySave(tmpMatrix.path, matrix)
  }
  return(tmpMatrix.path)
}

#' Title
#'
#' @param nmf.exp
#' @param k.min
#' @param k.max
#' @param outer.iter
#' @param inner.iter
#' @param conver.test.niter
#' @param conver.test.stop.threshold
#' @param out.dir
#'
#' @return
#'
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom data.table fread
#' @export
#'
#' @examples
runNmfGpu <- function(nmf.exp, k.min= 2, k.max = 2, outer.iter = 10,
                      inner.iter = 10^4, conver.test.niter = 10,
                      conver.test.stop.threshold = 40, out.dir = NULL,
                      tmp.path = "/tmp/nmf_tmp") {
  # Write raw matrix to tmp file
  tmpMatrix.path <- writeTmpMatrix(assay(nmf.exp, "raw"), tmp.path = tmp.path)
  # Define pattern to finde GPU_NMF output.
  tmp.dir <- dirname(tmpMatrix.path)
  h.pattern <- sprintf("%s_H.txt", basename(tmpMatrix.path))
  w.pattern <- sprintf("%s_W.txt", basename(tmpMatrix.path))
  # Check if deconv. matrix should be saved
  if(!is.null(out.dir)) dir.create(out.dir)
  # RUN NMF.
  dec.matrix <- lapply(k.min:k.max, function(k) {
    print(Sys.time())
    cat("Factorization rank: ", k, "\n")
    k.matrix <- lapply(1:outer.iter, function(i) {
      if(i %% 10 == 0) { cat("\tIteration: ", i, "\n") }
      frob.error <- 1
      while(frob.error == 1 | is.na(frob.error)){
        # VERSION 1.0
        # nmf.cmd <-
        #   sprintf('NMF_GPU %s -k %s -i %s', tmpMatrix.path, k, inner.iter)
        # nmf.stdout <- system(nmf.cmd, intern = T)
        nmf.cmd <- sprintf("%s -k %i -i %i -j %i -t %i", tmpMatrix.path, k,
                           inner.iter, conver.test.niter,
                           conver.test.stop.threshold)
        nmf.stdout <-
          system2("NMF_GPU", args = nmf.cmd, stdout = T, stderr = NULL)
        frob.error <- nmf.stdout[grep(nmf.stdout, pattern = "Distance")]
        frob.error <- as.numeric(sub(".*: ", "", frob.error))
        if(length(frob.error) == 0) frob.error <- 1
      }
      h.file <- list.files(tmp.dir, pattern = h.pattern, full.names = T)
      h.matrix <- fread(h.file, header = F)
      w.file <- list.files(tmp.dir, pattern = w.pattern, full.names = T)
      w.matrix <- fread(w.file, header = F)
      if(!is.null(out.dir)) {
        file.copy(h.file,
                  file.path(out.dir, sprintf("H_k%s_iter%s.txt", k, i)))
        file.copy(w.file,
                  file.path(out.dir, sprintf("W_k%s_iter%s.txt", k, i)))
      }
      file.remove(c(h.file, w.file))
      return(list(H = h.matrix,
                  W = w.matrix,
                  Frob.error = frob.error))
    })
    names(k.matrix) <- 1:outer.iter
    return(k.matrix)
  })
  names(dec.matrix) <- k.min:k.max
  ### Add NMF results to summarizedExp object.
  # Frob Errors
  frob.errors <- DataFrame(getFrobError(dec.matrix))
  colnames(frob.errors) <- as.character(k.min:k.max)
  nmf.exp <- setFrobError(nmf.exp, frob.errors)
  # H-Matrix List
  nmf.exp <- setHMatrixList(nmf.exp, getHMatrixList(dec.matrix))
  # W-Matrix List
  nmf.exp <- setWMatrixList(nmf.exp, getWMatrixList(dec.matrix))
  return(nmf.exp)
}

#' Title
#'
#' @param nmf.exp
#' @param k.min
#' @param k.max
#' @param outer.iter
#' @param inner.iter
#' @param conver.test.niter
#' @param conver.test.stop.threshold
#' @param out.dir
#'
#' @return
#'
#' @import SummarizedExperiment
#' @import NMF
#' @importFrom data.table fread
#' @importFrom S4Vectors DataFrame
#' @importFrom RcppCNPy npyLoad
#' @export
#'
#' @examples
runNmfGpuPyCuda <- function(nmf.exp, k.min= 2, k.max = 2, outer.iter = 10,
                            inner.iter = 10^4, conver.test.niter = 10,
                            conver.test.stop.threshold = 40, out.dir = NULL,
                            tmp.path = "/tmp/nmf_tmp", nmf.type = "N",
                            w.sparsness = 0, h.sparsness = 0, gpu.id = 0,
                            seed = FALSE, binary.file.transfer = FALSE,
                            cpu = FALSE) {
  if(!cpu){
    # Write raw matrix to tmp file
    tmpMatrix.path <- writeTmpMatrix(assay(nmf.exp, "raw"), tmp.path = tmp.path,
                                     binary = binary.file.transfer)
    # Define encoding
    encoding <- "txt"
    if(binary.file.transfer)
      encoding <- "npy"
    # Define pattern to finde GPU_NMF output.
    tmp.dir <- dirname(tmpMatrix.path)
    h.pattern <- sprintf("%s_H.%s", basename(tmpMatrix.path), encoding)
    w.pattern <- sprintf("%s_W.%s", basename(tmpMatrix.path), encoding)
    # Check if deconv. matrix should be saved
    if(!is.null(out.dir)) dir.create(out.dir)
    # RUN NMF.
    dec.matrix <- lapply(k.min:k.max, function(k) {
      print(Sys.time())
      cat("Factorization rank: ", k, "\n")
      k.matrix <- lapply(1:outer.iter, function(i) {
        if(i %% 10 == 0) { cat("\tIteration: ", i, "\n") }
        frob.error <- 1
        if (seed){
          nmf.cmd <-
            sprintf(paste0("%s %s -k %i -i %i -s %s -wo %s -ho %s ",
                           "-g %i -e %s -sets %s -sv %i"),
                    file.path(system.file(package = "Bratwurst"),
                              "python/nmf_mult.py"),
                    tmpMatrix.path, k, inner.iter, nmf.type, w.sparsness,
                    h.sparsness, gpu.id, encoding, "True", i)
          nmf.stdout <- system2("python", args = nmf.cmd, stdout = T,
                                stderr = NULL)
          frob.error <- nmf.stdout[grep(nmf.stdout, pattern = "Distance")]
          frob.error <- as.numeric(sub(".*: ", "", frob.error))
        } else {
          count <-  1
          while((frob.error == 1 | is.na(frob.error)) & count < 11){
            count <- count + 1
            nmf.cmd <-
              sprintf(paste0("%s %s -k %i -i %i -s %s -wo %s -ho %s ",
                             "-g %i -e %s"),
                      file.path(system.file(package = "Bratwurst"),
                                "python/nmf_mult.py"),
                      tmpMatrix.path, k, inner.iter, nmf.type, w.sparsness,
                      h.sparsness, gpu.id, encoding)
            nmf.stdout <- system2("python", args = nmf.cmd, stdout = T,
                                  stderr = NULL)
            frob.error <- nmf.stdout[grep(nmf.stdout, pattern = "Distance")]
            frob.error <- as.numeric(sub(".*: ", "", frob.error))
            if(length(frob.error) == 0) frob.error <- 1
          }
        }
        h.file <- list.files(tmp.dir, pattern = h.pattern, full.names = T)
        w.file <- list.files(tmp.dir, pattern = w.pattern, full.names = T)
        if(!binary.file.transfer){
          h.matrix <- fread(h.file, header = F)
          w.matrix <- fread(w.file, header = F)
        } else {
          h.matrix <- npyLoad(h.file)
          w.matrix <- npyLoad(w.file)
        }
        if(!is.null(out.dir)) {
          file.copy(h.file, file.path(out.dir, sprintf("H_k%s_iter%s.%s",
                                                       k, i, encoding)))
          file.copy(w.file, file.path(out.dir, sprintf("W_k%s_iter%s.%s",
                                                       k, i, encoding)))
        }
        file.remove(c(h.file, w.file))
        return(list(H = h.matrix,
                    W = w.matrix,
                    Frob.error = frob.error))
      })
      names(k.matrix) <- 1:outer.iter
      return(k.matrix)
    })
    names(dec.matrix) <- k.min:k.max
    ### Add NMF results to summarizedExp object.
    # Frob Errors
    frob.errors <- DataFrame(getFrobError(dec.matrix))
    colnames(frob.errors) <- as.character(k.min:k.max)
    nmf.exp <- setFrobError(nmf.exp, frob.errors)
    # H-Matrix List
    nmf.exp <- setHMatrixList(nmf.exp, getHMatrixList(dec.matrix))
    # W-Matrix List
    nmf.exp <- setWMatrixList(nmf.exp, getWMatrixList(dec.matrix))
  } else {
    # run NMF on CPU via CRAN R package
    # cpu.nmf.list <- runNmfCpu(nmf.exp, k.min = k.min, k.max = k.max,
    #                           outer.iter = outer.iter)
    cpu.nmf.list <- runNmfCpu(nmf.exp, k.min = k.min, k.max = k.max,
                              nrun = outer.iter, seed = seed)
    ### Add NMF results to summarizedExp object
    # compute Frob Erros and add them to nmf.exp
    frob.errors <- DataFrame(computeCpuFrobErrors(nmf.exp, cpu.nmf.list,
                                                  k.min, k.max))
    colnames(frob.errors) <- as.character(k.min:k.max)
    nmf.exp <- setFrobError(nmf.exp, frob.errors)
    # H-Matrix List
    nmf.exp <- setHMatrixList(nmf.exp, getCpuHMatrixList(cpu.nmf.list,
                                                         k.min, k.max))
    # W-Matrix List
    nmf.exp <- setWMatrixList(nmf.exp, getCpuWMatrixList(cpu.nmf.list,
                                                         k.min, k.max))
  }
 return(nmf.exp)
}

#' Title
#'
#' @param nmf.exp
#' @param k.min
#' @param k.max
#' @param outer.iter
#'
#' @return
#'
#' @import NMF
#'
#' @examples
#runNmfCpu <- function(nmf.exp, k.min = 2, k.max = 2, outer.iter = 10){
runNmfCpu <- function(nmf.exp, k.min = 2, k.max = 2, seed = FALSE, ...){
    ## create ExpressionSet
  eset <- new("ExpressionSet", exprs = assay(nmf.exp))

  # loop over all k and outer.iter iterations
  cpu.nmf.list <- lapply(k.min:k.max, function(k){
    print(Sys.time())
    cat("Factorization rank: ", k, "\n")
    if(seed){
      #cpu.nmf <- nmf(eset, rank = k, method = "brunet", nrun = outer.iter,
      cpu.nmf <- nmf(eset, rank = k, ...,
                     .options = list(keep.all = TRUE), seed = 12)
    }else {
      #cpu.nmf <- nmf(eset, rank = k, method = "brunet", nrun = outer.iter,
      cpu.nmf <- nmf(eset, rank = k, ...,
                     .options = list(keep.all = TRUE))
    }

    return(cpu.nmf)
  })

  return(cpu.nmf.list)
}

#' Title
#'
#' @param nmf.exp
#' @param k.min
#' @param k.max
#' @param outer.iter
#'
#' @return
#'
#' @import NMF
#'
#' @examples
computeCpuFrobErrors <- function(nmf.exp, cpu.nmf.list, k.min, k.max){
  frob.errors <- do.call(cbind, lapply(seq(length(k.min:k.max)), function(ki){
    # get all W-matrices
    W.list <- lapply(cpu.nmf.list[[ki]], basis)
    # get all H-matrices
    H.list <- lapply(cpu.nmf.list[[ki]], coef)
    fes <- sapply(seq(length(W.list)), function(oi){
      fe <- norm(assay(nmf.exp) - (W.list[[oi]] %*% H.list[[oi]]), type = "F") /
        norm(assay(nmf.exp), type = "F")
    })
    return(fes)
  }))

  return(frob.errors)
}


#==============================================================================#
#                               Getter FUNCTIONS                               #
#==============================================================================#
#' Getter function for FrobError from NMF-GPU list
#'
#' @param dec.matrix
#'
#' @return
#'
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#'
#' @examples
getFrobError <- function(dec.matrix) {
  frob.errorList <- lapply(dec.matrix, function(dec.m){
    lapply(dec.m, function(m) m$Frob.error)
  })
  frob.errorMatrix <- melt(frob.errorList)
  frob.errorMatrix <- dcast(frob.errorMatrix, L2 ~ L1, value.var = "value")
  frob.errorMatrix <-
    frob.errorMatrix[order(as.numeric(frob.errorMatrix$L2)), ,drop=FALSE]
  frob.errorMatrix <- frob.errorMatrix[, -1, drop=FALSE]
  frob.errorMatrix <-
    frob.errorMatrix[, order(as.numeric(colnames(frob.errorMatrix))), drop=FALSE]
  frob.errorMatrix[frob.errorMatrix == 1] <- NA
  return(frob.errorMatrix)
}

#' Getter function for H-Matrix list from NMF-GPU list
#'
#' @param dec.matrix
#'
#' @return
#'
#' @examples
getHMatrixList <- function(dec.matrix) {
  # Extract all entries for H and return.
  H.list <- lapply(dec.matrix, function(k.matrix) {
    H <- lapply(k.matrix, function(m) as.matrix(m$H))
    return(H)
  })
  return(H.list)
}

#' Getter function for W-Matrix list from NMF-GPU list
#'
#' @param dec.matrix
#'
#' @return
#'
#' @examples
getWMatrixList <- function(dec.matrix) {
  # Extract all entries for W and return.
  W.list <- lapply(dec.matrix, function(k.matrix) {
    W <- lapply(k.matrix, function(m) as.matrix(m$W))
    return(W)
  })
  return(W.list)
}

#' Getter function for W-Matrix list from NMF-CPU list
#'
#' @param dec.matrix
#'
#' @return
#'
#' @import NMF
#'
#' @examples
getCpuWMatrixList <- function(cpu.nmf.list, k.min, k.max) {
  # Extract all entries for W and return.
  W.list <- lapply(cpu.nmf.list, function(cpu.nmf) {
    W <- lapply(cpu.nmf, basis)
    return(W)
  })
  names(W.list) <- as.character(k.min:k.max)
  return(W.list)
}

#' Getter function for H-Matrix list from NMF-CPU list
#'
#' @param dec.matrix
#'
#' @return
#'
#' @examples
getCpuHMatrixList <- function(cpu.nmf.list, k.min, k.max) {
  # Extract all entries for H and return.
  H.list <- lapply(cpu.nmf.list, function(cpu.nmf) {
    H <- lapply(cpu.nmf, coef)
    return(H)
  })
  names(H.list) <- as.character(k.min:k.max)
  return(H.list)
}

#==============================================================================#
#             Criteria for optimal factorization rank - FUNCTIONS              #
#==============================================================================#
#' Compute basic statistics for Frobenius Errors
#'
#' @param nmf.exp
#'
#' @return
#' @export
#'
#' @examples
computeFrobErrorStats <- function(nmf.exp) {
  frob.errorMatrix <- as.matrix(FrobError(nmf.exp))
  min.frobError <- apply(frob.errorMatrix, 2, function(x) min(x, na.rm = T))
  sd.frobError <- apply(frob.errorMatrix, 2, function(x) sd(x, na.rm = T))
  mean.frobError <- colMeans(frob.errorMatrix, na.rm = T)
  cv.frobError <- sd.frobError / mean.frobError
  frobError.data <- DataFrame("k" = as.numeric(names(min.frobError)),
                              "min" = min.frobError,
                              "mean" = mean.frobError,
                              "sd" = sd.frobError,
                              "cv" = cv.frobError)
  nmf.exp <- setOptKStats(nmf.exp, frobError.data)
  return(nmf.exp)
}

#' Compute p-value with t-test for running K
#'
#' @param frobError.matrix
#'
#' @return
#' @export
#'
#' @examples
computePValue4FrobError <- function(frobError.matrix) {
  n.col <- ncol(frobError.matrix)
  p <- mapply(1:(n.col - 1), 2:n.col, FUN = function(i, j) {
    p <- t.test(frobError.matrix[, i], frobError.matrix[, j])
    return(p$p.value)
  })
  #p <- p.adjust(unlist(p), method = 'BH')
  return(-log10(unlist(p)))
}

#' Cosine similarity
#'
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
cosineSim <- function(a, b){
  return(sum(a * b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ))
}

#' Cosine distance
#'
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
cosineDist <- function(a, b){
  return(1 - cosineSim(a, b))
}

#' Create distance matrix with cosine similarity
#'
#' @param in.matrix
#' @param in.dimension
#'
#' @return
#' @export
#'
#' @examples
cosineDiss <- function(in.matrix, in.dimension=2){
  if(in.dimension == 1) in.matrix <- t(in.matrix)
  cosineDist.list <- lapply(1:ncol(in.matrix), function(i.outCol) {
    cosine.dists <- unlist(lapply(1:ncol(in.matrix), function(i.innCol) {
      cosineDist(in.matrix[, i.outCol], in.matrix[, i.innCol])
    }))
  })
  diss.matrix <- do.call(rbind, cosineDist.list)
  return(round(diss.matrix, digits = 14))
}

#' Create distance matrix with cosine similarity with matrix operations
#'
#' @param in.matrix
#' @param in.dimension
#'
#' @return
#' @export
#'
#' @examples
cosineDissMat <- function(in.matrix, in.dimension=2){
  if(in.dimension == 1) in.matrix <- t(in.matrix)
  squaredVectorSum <- apply(in.matrix, 2, function(m) { sqrt(sum(m * m)) })
  squaredVectorProduct <- squaredVectorSum %*% t(squaredVectorSum)
  squaredInputSum <- t(in.matrix) %*% in.matrix
  # sum(a*b) for any a,b in M
  diss.matrix <- 1 - squaredInputSum / squaredVectorProduct
  # CosineDistance = 1 - CosineSimilarity
  return(round(diss.matrix, digits = 14))
}

#' Compute Alexandrov Criterion --> Silhoutte Width
#'
#' @param nmf.exp
#'
#' @return
#'
#' @importFrom cluster pam
#' @export
#'
#' @examples
computeSilhoutteWidth <- function(nmf.exp) {
  sil.vec <- lapply(WMatrixList(nmf.exp), function(WMatrix.list) {
    nan.bool <- sapply(WMatrix.list, function(m) !any(is.nan(m)))
    concat.matrix <- do.call(cbind, WMatrix.list[which(nan.bool)])
    dist.matrix <- cosineDissMat(as.matrix(concat.matrix))
    my.pam <- pam(dist.matrix, k = ncol(WMatrix.list[[1]]),  diss = T)
    sil.sum <- sum(my.pam$silinfo$widths[, "sil_width"])
    sil.mean <- mean(my.pam$silinfo$widths[, "sil_width"])
    return(DataFrame(sumSilWidth = sil.sum,
                     meanSilWidth = sil.mean))
  })
  sil.vec <- do.call(rbind, sil.vec)
  nmf.exp <- setOptKStats(nmf.exp, cbind(OptKStats(nmf.exp), sil.vec))
  return(nmf.exp)
}

#' Compute Cophenetic correlation coefficient, TO BE IMPROVED
#'
#' @param nmf.exp
#'
#' @return
#' @export
#'
#' @examples
computeCopheneticCoeff <- function(nmf.exp) {
  coph.coeff <- lapply(WMatrixList(nmf.exp), function(WMatrix.list) {
    nan.bool <- sapply(WMatrix.list, function(m) !any(is.nan(m)))
    concat.matrix <- do.call(cbind, WMatrix.list[which(nan.bool)])
    #concat.matrix <- do.call(cbind, WMatrix.list)
    dist.matrix <- cosineDissMat(as.matrix(concat.matrix))
    my.hclust <- hclust(as.dist(dist.matrix))
    dist.cophenetic <- as.matrix(cophenetic(my.hclust))
    # take distance matrices without diagonal elements
    diag(dist.matrix) <- NA
    dist.matrix <- dist.matrix[which(!is.na(dist.matrix))]
    diag(dist.cophenetic) <- NA
    dist.cophenetic <- dist.cophenetic[which(!is.na(dist.cophenetic))]
    return(cor(cbind(dist.cophenetic, dist.matrix))[1, 2])
  })
  coph.coeff <- DataFrame(copheneticCoeff = unlist(coph.coeff))
  nmf.exp <- setOptKStats(nmf.exp, cbind(OptKStats(nmf.exp), coph.coeff))
  return(nmf.exp)
}

#' Compute amari type distance between two matrices
#'
#' @param matrix.A,matrix.B of the same dimensionality
#'
#' @return The amari type distance of matrix.A & matrix.B according
#'        to [Wu et. al, PNAS 2016]
#'
#' @references \url{http://www.pnas.org/content/113/16/4290.long}
#'
#' @export
#'
#' @examples
amariDistance <- function(matrix.A, matrix.B) {
  K <- dim(matrix.A)[2]
  C <- cor(matrix.A, matrix.B)
  return(1 - (sum(apply(C, 1, max)) + sum(apply(C, 2, max))) / (2 * K))
}

#' Compute Amari Distances from [Wu et. al, PNAS 2016]
#'
#' @param nmf.exp
#'
#' @return The average Amari-type error for each k
#'
#' @references \url{http://www.pnas.org/content/113/16/4290.long}
#'
#' @export
#'
#' @examples
computeAmariDistances <- function(nmf.exp){
  distance.averages <- lapply(WMatrixList(nmf.exp), function(matrices) {
    B <- length(matrices)
    distances.list <- unlist(lapply(1:(B - 1), function(b) {
      distances <- lapply((b + 1):B, function(b.hat) {
        amariDistance(matrices[[b]], matrices[[b.hat]])
      })
    }))
    return(mean(distances.list[!is.na(distances.list)]))
    # is.na to exclude corrupted matrices
  })
  distance.averages <- DataFrame(meanAmariDist = unlist(distance.averages))
  nmf.exp <- setOptKStats(nmf.exp, cbind(OptKStats(nmf.exp), distance.averages))
  return(nmf.exp)
}

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

#==============================================================================#
#                         H-MATRIX ANALYSIS FUNCTIONS                          #
#==============================================================================#
# Compute unbaised signature names by applying a row k-means
# and classify signature according to cluster and highest cluster mean.
#' Compute unsupervised signature names from colData and H-Matrix for given K
#'
#' @param nmf.exp
#' @param k.opt
#'
#' @return
#' @export
#'
#' @examples
getSignatureNames <- function(nmf.exp, k.opt) {
  H <- HMatrix(nmf.exp, k = k.opt)
  sig.names <- lapply(1:nrow(H), function(i) {
    k.cluster <- kmeans(as.numeric(H[i, ]), 2)
    n.kCluster <- which(k.cluster$centers == max(k.cluster$centers))
    samples.kCluster <- which(k.cluster$cluster == n.kCluster)
    meta.data <- colData(nmf.exp)[samples.kCluster, ]
    n.col <- ncol(meta.data)
    if (n.col > 1) {
      sigName.combs <- apply(as.data.frame(meta.data), 1, function(x){
        paste(x[2:n.col], collapse = " ")
      })
      sigName.combs <- sort(table(sigName.combs), decreasing = T)
      sig.name <- names(sigName.combs)[1]
      # TODO: Might be improvable by using exposures from H-Matrix.
      # Exposure proportion computation!?
      if(length(sigName.combs) > 1) {
        sig.prop <- round(sigName.combs / sum(sigName.combs), 2)
        sig.prop <- paste(names(sig.prop), sig.prop)
        sig.name <- paste(sig.prop, collapse = "\n")
      }
    } else {
      sig.name <- sprintf("Signature %s", i)
    }
    return(sig.name)
  })
  return(unlist(sig.names))
}


#==============================================================================#
#                         W-MATRIX ANALYSIS FUNCTIONS:                         #
#                             FEATURE SELECTION                                #
#==============================================================================#
#' Compute 'shannon' entropy per region.
#' High Entropy means highly specific for one signature.
#'
#' @param matrix
#'
#' @return
#'
#' @examples
computeEntropy <- function(matrix) {
  matrix.relativ <- t(apply(matrix, 1, function(x) x / sum(x)))
  matrix.entropy <- apply(matrix.relativ, 1, function(x) {
    p <- x * log2(length(x) * x)
    p[is.na(p)] <- 0
    h <- sum(p)
    return(h)
  })
  return(matrix.entropy)
}

#' Computes entropy for each feature in optimal K W-matrix
#'
#' @param nmf.exp
#'
#' @return
#' @export
#'
#' @examples
computeEntropy4OptK <- function(nmf.exp) {
  W <- WMatrix(nmf.exp, OptK(nmf.exp))
  W.entropy <- computeEntropy(W)
  nmf.exp <- setFeatureStats(nmf.exp, DataFrame("entropy" = W.entropy))
  return(nmf.exp)
}

#' Compute delta between each column (signature) per row
#'
#' @param matrix
#'
#' @return
#'
#' @examples
computeAbsDelta <- function(matrix) {
  delta.regions <- lapply(1:ncol(matrix), function(k) {
    delta.vec <- matrix[, k] - (rowSums(matrix[, -k]))
    return(delta.vec)
  })
  delta.regions <- do.call(cbind, delta.regions)
  return(delta.regions)
}

#' Computes absolut delta per feature for each signature given optimal K
#'
#' @param nmf.exp
#'
#' @return
#' @export
#'
#' @examples
computeAbsDelta4OptK <- function(nmf.exp) {
  W <- WMatrix(nmf.exp, OptK(nmf.exp))
  W.absDelta <- computeAbsDelta(W)
  nmf.exp <- setFeatureStats(nmf.exp, DataFrame("absDelta" = W.absDelta))
  return(nmf.exp)
}

#' Compute the coeficient of variation per row in a matrix.
#'
#' @param matrix
#'
#' @return
#'
#' @examples
computeCoefVar <- function(matrix) {
  apply(matrix, 1, function(r) sd(r) / mean(r))
}

#' Perform Kmeans on rows of a matrix to classify them into column (singature) combinations
#'
#' @param matrix
#'
#' @return
#'
#' @importFrom cluster silhouette
#' @export
#'
#' @examples
performRowKmeans <- function(matrix) {
  k.row <- apply(matrix, 1, function(x) {
    # Perform sigmoidal transformation to achieve better clustering
    x.trans <- sigmoidTransform(x)
    k.cluster <- kmeans(x.trans, 2)
    d <- dist(as.data.frame(x.trans))
    sil.mean <- mean(silhouette(k.cluster$cluster, d)[, 3])
    cluster.deltaMean <- mean(x[k.cluster$cluster == 1]) -
      mean(x[k.cluster$cluster == 2])
    return(list(centers = t(k.cluster$centers),
                silmean = sil.mean,
                explainedVar = k.cluster$betweenss / k.cluster$totss,
                oddsVar = sum(k.cluster$withinss) / k.cluster$betweenss,
                attribution = k.cluster$cluster,
                deltaMean = cluster.deltaMean))
  })
  return(k.row)
}

#' Compute Signature Combinations for W-matrix given optimal K
#'
#' @param nmf.exp
#' @param var.thres
#'
#' @return
#' @export
#'
#' @examples
computeFeatureStats <- function(nmf.exp, var.thres = 0.25) {
  W <- WMatrix(nmf.exp, OptK(nmf.exp))
  # Determine Region contributions
  k.row <- performRowKmeans(W)
  # Get kmeans stats
  k.explainedVar <- unlist(lapply(k.row, function(r) r$explainedVar))
  k.oddsVar <- unlist(lapply(k.row, function(r) r$oddsVar))
  # Extract mean silhouette width for each row.
  k.silmean <- unlist(lapply(k.row, function(r) r$silmean))
  # Compute delta between cluster centers (transformed data!)
  k.deltaCenter <-
    unlist(lapply(k.row, function(r) r$centers[1] - r$centers[2]))
  # Extract difference of cluster means
  k.deltaMean <- unlist(lapply(k.row, function(r) r$deltaMean))
  # Extract Signature combinations generated by k-means
  k.attribution <- lapply(k.row, function(r) abs(r$attribution))
  k.attribution <- do.call(rbind, k.attribution)
  k.ids <- apply(k.attribution, 1, function(r) paste(r, collapse = ""))
  # Determine regions which contribute to all signatures.
  k.varCoef <- computeCoefVar(W)
  all.signature <- which(k.varCoef < var.thres)
  k.ids[all.signature] <- gsub("2", "1", k.ids[all.signature])
  # Set FeatureStats
  feature.stats <- DataFrame("cluster" = k.ids,
                             "deltaCenters" = k.deltaCenter,
                             "deltaMean" = k.deltaMean,
                             "explainedVar" = k.explainedVar,
                             "oddsVar" = k.oddsVar,
                             "coefVar" = k.varCoef,
                             "meanSil" = k.silmean)
  # Map reverse cluster definitions to each other,
  # including sign change for delta means
  ids <- sort(unique(k.ids))[-1] ## this is not robust!!
  ## must be replaced by
  # equalPattern <- paste(rep("1", ncol(W)), collapse = "")
  # ids <- setdiff(sort(unique(k.ids)), equalPattern)
  ids1 <- ids[1:(length(ids) / 2)]
  ids2 <-  gsub("0", "2", gsub("2", "1", gsub("1", "0", ids1)))
  conv.id <- data.frame("id1" = ids1, "id2" = ids2)
  i.conv <- match(feature.stats$cluster, conv.id$id2)
  feature.stats$cluster[!is.na(i.conv)] <-
    as.character(conv.id$id1[i.conv[!is.na(i.conv)]])
  feature.stats$deltaMean[!is.na(i.conv)] <-
    (-1) * feature.stats$deltaMean[!is.na(i.conv)]
  feature.stats$deltaCenters[!is.na(i.conv)] <-
    (-1) * feature.stats$deltaCenters[!is.na(i.conv)]
  # Re-write cluster ids in a more useful binary code
  i <- which(feature.stats$deltaCenters > 0 &
               grepl(feature.stats$cluster, pattern = "2"))
  feature.stats$cluster[i] <- gsub("2", "0", feature.stats$cluster[i])
  i <- which(feature.stats$deltaCenters < 0 &
               grepl(feature.stats$cluster, pattern = "2"))
  feature.stats$cluster[i] <-
    gsub("2", "1", gsub("1", "0", feature.stats$cluster[i]))
  # Add featureStats to NMF experiment
  nmf.exp <- setFeatureStats(nmf.exp, feature.stats)
  return(nmf.exp)
}

#' Title
#'
#' @param k.ids
#'
#' @return
#' @export
#'
#' @examples
getIndex4SignatureCombs <- function(k.ids) {
  # Get for all signature combinations corresponding regions.
  ids <- sort(unique(k.ids))[-1]
  n.ids <- length(ids) + 1
  i.regions <- lapply(1:(n.ids / 2), function(i) {
    sub.id <- ids[c(i, n.ids - i)]
    i.region <- which(k.ids %in% sub.id)
    return(i.region)
  })
  names(i.regions) <- ids[1:(n.ids / 2)]
  # Get regions for all signature id.
  allSig.id <- sort(unique(k.ids))[1]
  i.regions[[allSig.id]] <- which(k.ids %in% allSig.id)
  return(i.regions)
}

#' Title
#'
#' @param W
#' @param i.regions
#'
#' @return
#' @export
#'
#' @examples
computeMeanDiff4SignatureCombs <- function(W, i.regions) {
  # Compute Row mean diff between Cluster 1 and 2.
  w.diffs <- lapply(1:length(i.regions), function(i) {
    w <- W[i.regions[[i]], ]
    if (!grepl(names(i.regions)[i], pattern = "2")) {
      w.diff <- rowMeans(w)
    } else {
      signature.comb <-
        as.numeric(unlist(strsplit(names(i.regions)[[i]], split = "")))
      w.mean1 <- rowMeans(w[, which(signature.comb == 1)])
      w.mean2 <- rowMeans(w[, which(signature.comb == 2)])
      w.diff <- w.mean1 - w.mean2
    }
    return(w.diff)
  })
  names(w.diffs) <- names(i.regions)
}

### Get number of peaks per signature combination and create a barplot.
#' Title
#'
#' @param i.regions
#' @param w.diffs
#'
#' @return
#' @export
#'
#' @examples
getSignatureCombCounts <- function(i.regions, w.diffs) {
  n.peaks <- lapply(names(i.regions), function(region.n) {
    n.peak1 <- sum(w.diffs[[region.n]] > 0)
    n.peak2 <- sum(w.diffs[[region.n]] < 0)
    n.peak <- c(n.peak1, n.peak2)
    n.peak <- melt(n.peak)
    n.peak$clusterId <- c("1", "2")
    n.peak$sigCombId <- region.n
    return(n.peak)
  })
  n.peaks <- do.call(rbind, n.peaks)
  return(n.peaks)
}


#' Normalize the signatures matrix (W)
#'
#' After column normalization of the matrix W, the inverse factors are
#' mutiplied with the rows of H in order to keep the matrix product W*H
#' constant.
#'
#' @param nmf.exp
#'
#' @return A data structure of type nmfExperiment
#'
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
#' @examples
#'  NULL
#'
normalizeW <- function(nmf.exp){
  # account for WMatrixList and HMatrixList
  all_list <- lapply(seq_along(WMatrixList(nmf.exp)), function(k_ind){
    k_list <-
      lapply(seq_along(WMatrixList(nmf.exp)[[k_ind]]), function(init_ind){
        tempW <- WMatrixList(nmf.exp)[[k_ind]][[init_ind]]
        tempH <- HMatrixList(nmf.exp)[[k_ind]][[init_ind]]
        normFactor <- colSums(tempW)
        # catch errors associated with NaNs in W or H
        if (any(is.nan(normFactor))){
          return(list(W = tempW,
                      H = tempH))
        }else{
          newSigs <- as.matrix(normalize_df_per_dim(tempW, 2))
          newExpo <- tempH * normFactor
          #newV <- newSigs %*% newExpo
          #oldV <- tempW %*% tempH
          return(list(W = newSigs,
                      H = newExpo))
        }
      })
    names(k_list) <- names(WMatrixList(nmf.exp)[[k_ind]])
    return(k_list)
  })
  names(all_list) <- names(WMatrixList(nmf.exp))
  thisWMatrixList <- lapply(all_list, function(current_k_list){
    kWMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$W)
    })
  })
  nmf.exp <- setWMatrixList(nmf.exp, thisWMatrixList)
  thisHMatrixList <- lapply(all_list, function(current_k_list){
    kHMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$H)
    })
  })
  nmf.exp <- setHMatrixList(nmf.exp, thisHMatrixList)
  return(nmf.exp)
}

#' Normalize the signatures matrix (H)
#'
#' After row normalization of the matrix H, the inverse factors are
#' mutiplied with the columns of W in order to keep the matrix product W*H
#' constant.
#'
#' @param nmf.exp
#'
#' @return A data structure of type nmfExperiment
#'
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
#' @examples
#'  NULL
#'
normalizeH <- function(nmf.exp){
  # account for WMatrixList and HMatrixList
  all_list <- lapply(seq_along(WMatrixList(nmf.exp)), function(k_ind){
    k_list <-
      lapply(seq_along(WMatrixList(nmf.exp)[[k_ind]]), function(init_ind){
        tempW <- WMatrixList(nmf.exp)[[k_ind]][[init_ind]]
        tempH <- HMatrixList(nmf.exp)[[k_ind]][[init_ind]]
        normFactor <- rowSums(tempH)
        newExpo <- as.matrix(normalize_df_per_dim(tempH, 1))
        newSigs <- tempW * normFactor
        return(list(W = newSigs,
                    H = newExpo))
      })
    names(k_list) <- names(WMatrixList(nmf.exp)[[k_ind]])
    return(k_list)
  })
  names(all_list) <- names(WMatrixList(nmf.exp))
  thisWMatrixList <- lapply(all_list, function(current_k_list){
    kWMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$W)
    })
  })
  nmf.exp <- setWMatrixList(nmf.exp, thisWMatrixList)
  thisHMatrixList <- lapply(all_list, function(current_k_list){
    kHMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$H)
    })
  })
  nmf.exp <- setHMatrixList(nmf.exp, thisHMatrixList)
  return(nmf.exp)
}

#' Regularize the signatures matrix (H)
#'
#' After row regularization of the matrix H, the inverse factors are
#' mutiplied with the columns of W in order to keep the matrix product W*H
#' constant.
#'
#' @param nmf.exp
#'
#' @return A data structure of type nmfExperiment
#'
#' @export
#'
#' @examples
#'  NULL
#'
regularizeH <- function(nmf.exp){
  # account for WMatrixList and HMatrixList
  all_list <- lapply(seq_along(WMatrixList(nmf.exp)), function(k_ind){
    k_list <-
      lapply(seq_along(WMatrixList(nmf.exp)[[k_ind]]), function(init_ind){
        tempW <- WMatrixList(nmf.exp)[[k_ind]][[init_ind]]
        tempH <- HMatrixList(nmf.exp)[[k_ind]][[init_ind]]
        normFactor <- rowMax(tempH)
        newExpo <- tempH / normFactor
        newSigs <- tempW * normFactor
        return(list(W = newSigs,
                    H = newExpo))
      })
    names(k_list) <- names(WMatrixList(nmf.exp)[[k_ind]])
    return(k_list)
  })
  names(all_list) <- names(WMatrixList(nmf.exp))
  thisWMatrixList <- lapply(all_list, function(current_k_list){
    kWMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$W)
    })
  })
  nmf.exp <- setWMatrixList(nmf.exp, thisWMatrixList)
  thisHMatrixList <- lapply(all_list, function(current_k_list){
    kHMatrixList <- lapply(current_k_list, function(current_entry){
      return(current_entry$H)
    })
  })
  nmf.exp <- setHMatrixList(nmf.exp, thisHMatrixList)
  return(nmf.exp)
}

#' Merge two objects of type nmfExperiment
#'
#' @param in_nmf1 First object of type nmfExperiment
#' @param in_nmf2 Second object of type nmfExperiment
#' @param rerunStats
#'  Boolean to indicate whether statistics should be computed on the merged
#'  object
#' @param verbose
#'
#' @return The merged object of type nmfExperiment
#' @export
#'
#' @examples
#'  NULL
#'
merge.nmf <- function(in_nmf1,
                      in_nmf2,
                      rerunStats = TRUE,
                      verbose = FALSE){
  error_msg <- "Bratwurst:::merge.nmf::error: input type mismatch.\n"
  ## check if the two stem from the same input
  if(all(names(in_nmf1) == names(in_nmf2)) &
     length(assays(in_nmf1)) == length(assays(in_nmf2))){
    if(verbose) cat("Bratwurst:::merge.nmf::verbose:Features of the two",
                    "input data structures match.\n")
    assayLogicVector <-
      unlist(lapply(seq_along(assays(in_nmf1)), function(current_ind){
        temp1 <- assays(in_nmf1)[[current_ind]]
        temp2 <- assays(in_nmf2)[[current_ind]]
        all(temp1 == temp2)
      }))
    if(all(assayLogicVector) &
       all(names(HMatrixList(in_nmf1)) == names(HMatrixList(in_nmf2))) &
       all(names(WMatrixList(in_nmf1)) == names(WMatrixList(in_nmf2)))){
      if(verbose) cat("Bratwurst:::merge.nmf::verbose:Assays of the two",
                      "input data structures match.\n")
      ## initialize
      if(verbose) cat("Bratwurst:::merge.nmf::verbose:Initialize.\n")
      merged.nmf.exp <- in_nmf1
      ## merge HMatrixList
      if(verbose) cat("Bratwurst:::merge.nmf::verbose:Concatenate",
                      "HMatrixList.\n")
      temp_list <-
        lapply(names(HMatrixList(in_nmf1)), function(currentRank){
          c(HMatrixList(in_nmf1, k = currentRank),
            HMatrixList(in_nmf2, k = currentRank))
        })
      names(temp_list) <- names(HMatrixList(in_nmf1))
      merged.nmf.exp <- setHMatrixList(merged.nmf.exp, temp_list)
      ## merge WMatrixList
      if(verbose) cat("Bratwurst:::merge.nmf::verbose:Concatenate",
                      "WMatrixList.\n")
      temp_list <-
        lapply(names(WMatrixList(in_nmf1)), function(currentRank){
          c(WMatrixList(in_nmf1, k = currentRank),
            WMatrixList(in_nmf2, k = currentRank))
        })
      names(temp_list) <- names(WMatrixList(in_nmf1))
      merged.nmf.exp <- setWMatrixList(merged.nmf.exp, temp_list)
      ## merge FrobError
      if(verbose) cat("Bratwurst:::merge.nmf::verbose:Concatenate",
                      "FrobError.\n")
      merged.nmf.exp@FrobError <- rbind(FrobError(in_nmf1), FrobError(in_nmf2))
      ## recalculate OptKStats if already available
      if(verbose) cat("Bratwurst:::merge.nmf::verbose:Recalculate error",
                      "statistics if necessary.\n")
      # if(sum(dim(OptKStats(in_nmf1))) + sum(dim(OptKStats(in_nmf2))) > 0){
      if(rerunStats){
        ## recalculate FrobErrorStats
        if(verbose) cat("Bratwurst:::merge.nmf::verbose:Recalculate",
                        "FrobErrorStats.\n")
        merged.nmf.exp <- computeFrobErrorStats(merged.nmf.exp)
        ## recalculate Alexandrov Criterion
        if(verbose) cat("Bratwurst:::merge.nmf::verbose:Recalculate",
                        "Alexandrov Criterion.\n")
        merged.nmf.exp <- computeSilhoutteWidth(merged.nmf.exp)
        # recalculate Cophenetic correlation coefficient
        if(verbose) cat("Bratwurst:::merge.nmf::verbose:Recalculate",
                        "Cophenetic correlation coefficient.\n")
        merged.nmf.exp <- computeCopheneticCoeff(merged.nmf.exp)
        # recalculate Amari type distance
        if(verbose) cat("Bratwurst:::merge.nmf::verbose:Recalculate",
                        "Amari type distance.\n")
        merged.nmf.exp <- computeAmariDistances(merged.nmf.exp)
      }
      return(merged.nmf.exp)
    } else {
      cat(error_msg)
      return(NULL)
    }
  } else {
    cat(error_msg)
    return(NULL)
  }
}

#' Compute signature specific features
#'
#' @param nmf.exp A NMF experiment object
#' @param rowDataId The index of the rowData(nmf.exp) data.frame that should be used for
#'  feature extraction. In case rowData(nmf.exp)[,1] is a GRanges or a related object like
#'  GenomicInteractions this parameter can be ignored
#'
#' @return nmf.exp with filles SignatureSpecificFeatures container
#'
#'
#' @export
#'
#' @examples
#'
computeSignatureSpecificFeatures <- function(nmf.exp, rowDataId = 3){
  if (length(OptK(nmf.exp)) == 0){
    stop("You need to first define an optimal k before being able to compute
           signature specific features!")
  } else {
    if (nrow(FeatureStats(nmf.exp)) == 0){
      message("Computing feature stats...")
      nmf.exp <- computeFeatureStats(nmf.exp)
    }
    fstats <- FeatureStats(nmf.exp)
    # identify unique cluster membership strings
    clusterMemberships <-
      sapply(unique(as.character(fstats$cluster)),
             function (x) lengths(regmatches(x, gregexpr("1", x))))
    sigSpecClusters <-
      sort(names(clusterMemberships[which(clusterMemberships == 1)]),
           decreasing = TRUE)

    if (class(rowData(nmf.exp)[, 1]) %in% c("GRanges", "GInteractions",
                                           "GenomicInteractions")){
      signatureSpecificFeatures <- lapply(sigSpecClusters, function(ssc){
        features <- rowData(nmf.exp)[, 1][which(fstats$cluster == ssc)]
        return(features)
      })
    }else{
      signatureSpecificFeatures <- lapply(sigSpecClusters, function(ssc){
        features <- rowData(nmf.exp)[, rowDataId][which(fstats$cluster == ssc)]
        return(features)
      })
    }
    names(signatureSpecificFeatures) <- sigSpecClusters
    nmf.exp@SignatureSpecificFeatures <- signatureSpecificFeatures
    return(nmf.exp)
  }
}


#==============================================================================#
#                      NMF tensorflow Wrapper - FUNCTION                       #
#==============================================================================#

#' Computes NMF on tensorflow using the reticulate framework
#'
#' @param nmf.exp
#' @param k.min
#' @param k.max
#' @param outer.iter
#' @param inner.iter
#' @param conver.test.niter
#' @param conver.test.stop.threshold
#'
#' @return
#'
#' @import SummarizedExperiment
#' @importFrom data.table fread
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
runNMFtensor <- function (nmf.exp, k.min = 2, k.max = 2, outer.iter = 10, inner.iter = 10^4, 
                          conver.test.niter = 10, conver.test.stop.threshold = 40){
  # Convert params to integer
  nmf.params <- lapply(list(k.min=k.min, 
                            k.max=k.max, 
                            outer.iter=outer.iter, 
                            inner.iter=inner.iter, 
                            conver.test.niter=conver.test.niter, 
                            conver.test.stop.threshold=conver.test.stop.threshold), 
                       as.integer)
  
  # Source NMF tensorflow python script
  source_python(file.path(system.file(package = "Bratwurst"), "python/nmf_tensor.py"))
  
  # Run NMF
  X <- assay(nmf.exp, "raw")
  dec.matrix <- lapply(k.min:k.max, function(k) {
    k <- as.integer(k)
    
    print(Sys.time())
    cat("Factorization rank: ", k, "\n")
    k.matrix <- lapply(1:outer.iter, function(i) {
      if (i%%10 == 0) cat("\tIteration: ", i, "\n")
      
      nmf.eval <- NMF(X, k, nmf.params$inner.iter, nmf.params$conver.test.niter, nmf.params$conver.test.stop.threshold)
      names(nmf.eval) <- c("W", "H", "iterations", "Frob.error")
      return(nmf.eval)
    })
    names(k.matrix) <- 1:outer.iter
    
    iters <- paste(sapply(k.matrix, function(x) {x$iterations}), collapse = ",")
    print(paste("NMF converged after ", iters, "iterations"))
    
    return(k.matrix)
  })
  
  # Build NMF object
  names(dec.matrix) <- k.min:k.max
  frob.errors <- DataFrame(getFrobError(dec.matrix))
  colnames(frob.errors) <- as.character(k.min:k.max)
  nmf.exp <- setFrobError(nmf.exp, frob.errors)
  nmf.exp <- setHMatrixList(nmf.exp, getHMatrixList(dec.matrix))
  nmf.exp <- setWMatrixList(nmf.exp, getWMatrixList(dec.matrix))
  return(nmf.exp)
}

