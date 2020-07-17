#' Microarray expression data of B-ALL, T-ALL and AML samples.
#'
#' Dataset containing Affymetrix Hgu6800 microarray expression data of B-ALL,
#' T-ALL and AML samples.
#' This dataset had also been used by Brunet et al. (PNAS, 2004) and
#' Gaujoux et al. (BMC Bioinformatics, 2010).
#'
#' @format A list with an expression matrix and annotatin data frame
#' \describe{
#'   \item{matrix}{4452 rows/genes and 38 samples}
#'   \item{annotation}{AL and _AML status with type}
#'   ...
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/10521349/}
#'
#' @usage data(leukemia)
#'
"leukemia"

#'
#'
#' #' Compute NMF on the test dataset
#' #'
#' #'
#' #' @return NMF object
#' #'
#' #'
#' #' @examples
#' #' \dontrun{
#' #' NMF_leukemia()
#' #' }
#' NMF_leukemia <- function() {
#'   data("leukemia")
#'   leukemia
#'   runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#'                     method = "NMF",
#'                     n_initializations = 10,
#'                     extract_features = TRUE)
#' }
