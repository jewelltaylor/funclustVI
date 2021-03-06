% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{get_approx_probability_matrix_hcl}
\alias{get_approx_probability_matrix_hcl}
\title{Generates probability matrix with approximated probability values for each cluster}
\usage{
get_approx_probability_matrix_hcl(Y, K, x)
}
\arguments{
\item{Y}{The matrix containing rows corresponding the curves}

\item{K}{The number of clusters in cluster data}

\item{x}{The x used to generate the clusters}
}
\value{
A probability matrix with probability values for each cluster
}
\description{
Generates probability matrix with approximated probability values for each cluster
}
\examples{
get_equal_probability_matrix(K, curves_per_cluster, x)
}
