# Copyright © 2015-2020  The ButchR package contributors
# This file is part of the ButchR package. The ButchR package is licenced
# under GPL-3

#' Function to normalizeUpperQuartile on a matrix
#'
#' @param matrix input matrix
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' h <- HMatrix(nmf_exp, k = 5)
#' normalizeUpperQuartile(h)
#' }
normalizeUpperQuartile <- function(matrix) {
  matrix.norm <- apply(matrix, 2, function(c) {
    nf <- c[c != 0]
    c <- c / quantile(nf, 0.75)
    return(c)
  })
  return(matrix.norm)
}

#' Function to used to perform a k-means k=2 on a W matrix
#'
#' @param col.vector vector taken from a matrix columns
#' @param q quantile
#'
#' @return
#'
#' @examples
#' \dontrun{
#' Wf <- HMatrix(nmf_exp, k = 5)
#'ssf <- apply(Wf, 1, function(x){
#' x <- sigmoidTransform(x)
#' k <- kmeans(x, 2)
#' max_idx <- which.max(k$centers)
#' paste(if_else(k$cluster == max_idx, "1", "0"), collapse = "")
#' })
#' }
sigmoidTransform <- function(col.vector, q = 0.95) {
  q <- as.numeric(quantile(col.vector, q))
  x <- 2 / (1 + exp((-2) * col.vector / q)) - 1
  return(x)
}

#' Function to used to perform a k-means k=2 on a W matrix - vversion 2
#'
#' @param col.vector vector taken from a matrix columns
#' @param q quantile
#'
#' @return
#'
#' @examples
#' \dontrun{
#' Wf <- HMatrix(nmf_exp, k = 5)
#'ssf <- apply(Wf, 1, function(x){
#' x <- sigmoidTransform2(x)
#' k <- kmeans(x, 2)
#' max_idx <- which.max(k$centers)
#' paste(if_else(k$cluster == max_idx, "1", "0"), collapse = "")
#' })
#' }
sigmoidTransform2 <- function(col.vector, q = 0.95) {
  q <- as.numeric(quantile(col.vector, q))
  x <- 1 / (1 + exp((-2) * col.vector / q))
  return(x)
}

#' Function to rankTransforma matrix
#'
#' @param matrix input matrix
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' h <- HMatrix(nmf_exp, k = 5)
#' rankTransform(h)
#' }
rankTransform  <- function(matrix) {
  trans.matrix <- apply(matrix, 2, function(c) {
    rank(c) / length(c)
  })
  return(trans.matrix)
}


#' Function to order binary data matrix
#'
#' @param matrix input matrix
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' orderBinary(bin_matrix)
#' }
orderBinary <- function(matrix) {
  col.sum <- apply(matrix, 2, sum)
  unlist(sapply(unique(col.sum), function(s) which(col.sum == s)))
}
