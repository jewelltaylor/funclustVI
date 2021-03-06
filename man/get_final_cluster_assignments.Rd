% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{get_final_cluster_assignments}
\alias{get_final_cluster_assignments}
\title{Gets a vector of the final cluster assignments based on the probability matrix}
\usage{
get_final_cluster_assignments(probability_matrix)
}
\arguments{
\item{probability_matrix}{A matrix in which the rows represent the probabilities that the curve is in each of the clusters}
}
\value{
A vector of the final cluster assignments
}
\description{
Gets a vector of the final cluster assignments based on the probability matrix
}
\examples{
get_final_cluster_assignments(probability_matrix)
}
