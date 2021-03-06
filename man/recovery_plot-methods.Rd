% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/plot_functions_lite.R
\docType{methods}
\name{recovery_plot}
\alias{recovery_plot}
\alias{recovery_plot,matrix-method}
\alias{recovery_plot,ANY,ANY-method}
\title{Recovery plot}
\usage{
recovery_plot(x, annot, ...)

\S4method{recovery_plot}{matrix}(x, annot)
}
\arguments{
\item{x}{an object of class matrix, ButchR_NMF, ButchR_joinNMF, or
ButchR_integrativeNMF.}

\item{annot}{input named factor/character with biological/clinical
annotation for the columns of the input data.}

\item{...}{additional parameters.}
}
\value{
recovery plot or a data frame with the association analysis results.
}
\description{
Generate a recovery plot to visualize the association between the NMF
signatures and a biological or clinical annotation factor.
}
\examples{
\dontrun{
recovery_plot(x, annot_factor)
}
\dontrun{
data(leukemia)
leukemia_nmf_exp <- run_NMF_tensor(X = leukemia$matrix,
                                   ranks = 2:4,
                                   method = "NMF",
                                   n_initializations = 10,
                                   extract_features = TRUE)
recovery_plot(HMatrix(leukemia_nmf_exp, k = 4),
              leukemia$annotation$ALL_AML)
}
}
