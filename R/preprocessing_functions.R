# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#' Title
#'
#' @param matrix
#'
#' @return
#' @export
#'
#' @examples
normalizeUpperQuartile <- function(matrix) {
  matrix.norm <- apply(matrix, 2, function(c) {
    nf <- c[c != 0]
    c <- c / quantile(nf, 0.75)
    return(c)
  })
  return(matrix.norm)
}

#' Title
#'
#' @param col.vector
#' @param q
#'
#' @return
#' @export
#'
#' @examples
sigmoidTransform <- function(col.vector, q = 0.95) {
  q <- as.numeric(quantile(col.vector, q))
  x <- 2 / (1 + exp((-2) * col.vector / q)) - 1
  return(x)
}

#' Title
#'
#' @param col.vector
#' @param q
#'
#' @return
#' @export
#'
#' @examples
sigmoidTransform2 <- function(col.vector, q = 0.95) {
  q <- as.numeric(quantile(col.vector, q))
  x <- 1 / (1 + exp((-2) * col.vector / q))
  return(x)
}

#' Title
#'
#' @param matrix
#'
#' @return
#' @export
#'
#' @examples
rankTransform  <- function(matrix) {
  trans.matrix <- apply(matrix, 2, function(c) {
    rank(c) / length(c)
  })
  return(trans.matrix)
}


#' Function to order binary data matrix
#'
#' @param matrix
#'
#' @return
#' @export
#'
#' @examples
orderBinary <- function(matrix) {
  col.sum <- apply(matrix, 2, sum)
  unlist(sapply(unique(col.sum), function(s) which(col.sum == s)))
}
