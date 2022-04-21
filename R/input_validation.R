#' Validate if input matrix for NMF
#'
#' @param X matrix
#' @return X
#' @examples
#' \dontrun{
#' ButchR:::val_nonnegative_matrix(X)
#' }
val_nonnegative_matrix <- function(X) {
  if(is.data.frame(X)) {
    warning("Provided input is a data frame, coercing into matrix.")
    X = as.matrix(X)
  }
  if (!is.matrix(X)) {
    stop("\nProvided input matrix is not a matrix\n")
  }
  if (!is.numeric(X[1,1])) {
    stop("\nProvided input matrix is not a numeric matrix\n")
  }
  if (min(X) < 0 ) {
    stop("\nNegative values present in input matrix\n
         only non-negative matrices supported\n")
  }
  return(X)
}


#' Validate input ranks
#'
#' @param ranks numeric with ranks
#' @param ncolmat ncol of input matrix
#' @return NULL
#' @examples
#' \dontrun{
#' ButchR:::val_ranks_torun(X)
#' }
val_ranks_torun <- function(ranks, ncolmat) {
  if (!is.numeric(ranks)) {
    stop("\nProvided factorization ranks is not a numeric vector\n")
  }
  if(max(ranks) >= ncolmat) {
    stop("\nMaximum factorization rank should be less than the number of columns
         in the input matrix\n")
  }
  if(min(ranks) < 2) {
    stop("\nMinimum factorization rank should be greater or equal than 2\n")
  }
  NULL
}

#' Validate input single numeric integer
#'
#' @param x numeric
#' @param id param id
#' @param minval minimum value to test
#' @return NULL
#' @examples
#' \dontrun{
#' ButchR:::val_single_integer(X)
#' }
val_single_integer <- function(x, id, minval = 0) {
  y <- paste(x, collapse=",")
  if (!is.numeric(x)) {
    stop("\n", id, " = ", y, ": expecting single numeric value\n")
  }
  if (!length(x) == 1) {
    stop("\n", id, " = ", y, ": expecting single numeric value\n")
  }
  if (!x%%1==0) {
    stop("\n", id, " = ", y, ": should be a positive integer\n")
  }
  if (x < minval) {
    stop("\n", id, " = ", y, ": should be an integer >= ", minval,"\n")
  }
  NULL
}


#' Validate input single numeric
#'
#' @param x numeric
#' @param id param id
#' @param minval minimum value to test
#' @return NULL
#' @examples
#' \dontrun{
#' ButchR:::val_single_numeric(X)
#' }
val_single_numeric <- function(x, id, minval = 0) {
  if (!is.numeric(x)) {
    stop("\n", id, " = ", x, ": expecting single numeric value\n")
  }
  if (!length(x) == 1) {
    stop("\n", id, " = ", x, ": expecting single numeric value\n")
  }
  if (x < minval ) {
    stop("\n", id, " = ", x, ": should be  >= ", minval,"\n")
  }
  NULL
}


#' Validate input positive numeric vector
#'
#' @param x numeric
#' @param id param id
#' @return NULL
#' @examples
#' \dontrun{
#' ButchR:::val_positive_numeric(X)
#' }
val_positive_numeric <- function(x, id) {
  if (!is.numeric(x)) {
    stop("\n", id, " = ", x, ": expecting a numeric vector\n")
  }
  if (length(x) < 1) {
    stop("\n", id, " = ", x, ": expecting a numeric vector of length >= 1\n")
  }
  if (!all(x >= 0 )) {
    stop("\n", id, " = ",
         paste0(x, collapse = ";"),
         ": should all be  >= 0\n")
  }
  NULL
}



#' Validate input single boolean
#'
#' @param x bool
#' @param id param id
#' @return NULL
#' @examples
#' \dontrun{
#' ButchR:::val_single_numeric(X)
#' }
val_single_boolean <- function(x, id) {
  if (!is.logical(x)) {
    stop("\n", id, " = ", x, ": expecting single boolean value: TRUE or FALSE\n")
  }
  if (!length(x) == 1) {
    stop("\n", id, " = ", x, ": expecting single boolean value: TRUE or FALSE\n")
  }
  NULL
}



#' Validate if input graph is a matrix for GRNMF_SC
#'
#' @param graph graph
#' @param X matrix
#' @param method NMF method
#' @return graph
#' @examples
#' \dontrun{
#' ButchR:::val_graph_GRNMF_SC(graph)
#' }
val_graph_GRNMF_SC <- function(graph, X, method) {
  if(!method == "GRNMF_SC") {
    warning("\ngraph: graph is only used with the method: GRNMF_SC.\n",
            "\nignoring graph\n")
  }
  if(is.data.frame(graph)) {
    warning("\ngraph: Provided input graph is a data frame, coercing into matrix.")
    graph = as.matrix(graph)
  }
  if (!is.matrix(graph)) {
    stop("\ngraph: Provided input graph is not a matrix\n")
  }
  if (!nrow(graph) == ncol(graph)) {
    stop("\ngraph: Expecting square matrix found: nrow(graph) != ncol(graph)!\n")
  }
  if (!ncol(X) == ncol(graph)) {
    stop("\ngraph: Expecting matrix with same number of columns than input X\n",
         "found: ncol(X)=", ncol(X), " and ncol(graph)=", ncol(graph), "\n")
  }
  if (!is.numeric(graph[1,1])) {
    stop("\ngraph: Provided input graph is not a numeric matrix\n")
  }
  # if (min(graph) < 0 ) {
  #   stop("\nNegative values present in input graph\n
  #        only non-negative graph matrices supported\n")
  # }
  return(graph)
}
