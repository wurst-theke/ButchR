% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\docType{methods}
\name{OptKStats}
\alias{OptKStats}
\title{Return optimal factorization rank (K) Statistics}
\usage{
OptKStats(x, ...)
}
\arguments{
\item{x}{an nmfExperiment_lite object}

\item{...}{additional parameters.}
}
\value{
optimal K Statistics
}
\description{
Return optimal factorization rank (K) Statistics
}
\examples{
\dontrun{
data("leukemia")
nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
method = "NMF",
n_initializations = 2)
OptKStats(nmf_exp)
}
}