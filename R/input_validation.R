#' Validate if input matrix for NMF
#'
#' @param X matrix
#' @return X
#' @examples
#' \dontrun{
#' ButchR:::val_nonnegative_matrix(X)
#' }
val_nonnegative_matrix <- function(X) {
  if (!is.numeric(X[1,1])) {
    stop("\nProvided input matrix is not a numeric matrix\n")
  }
  if(is.data.frame(X)) {
    warning("Provided input is a data frame, coercing into matrix.")
    X = as.matrix(X)
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
  if(min(ranks) <= 0) {
    stop("\nMinimum factorization rank should be greater or equal than 2\n")
  }
  NULL
}

#' Validate input single numeric
#'
#' @param x numeric
#' @param id param id
#' @return NULL
#' @examples
#' \dontrun{
#' ButchR:::val_single_numeric(X)
#' }
val_single_numeric <- function(x, id) {
  if (!is.numeric(x)) {
    stop("\n", id, " = ", x, ": expecting single numeric value\n")
  }
  if (!length(x) == 1) {
    stop("\n", id, " = ", x, ": expecting single numeric value\n")
  }
  if (x < 0 | !x%%1==0) {
    stop("\n", id, " = ", x, ": should be a positive integer\n")
  }
  NULL
}


