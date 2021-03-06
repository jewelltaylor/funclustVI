% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{update_m_list}
\alias{update_m_list}
\title{Update the m parameter for each cluster}
\usage{
update_m_list(
  Y,
  m_not_vector,
  phi_matrix,
  K,
  A_vector,
  R_vector,
  probability_matrix,
  x,
  nbasis,
  sigma_list
)
}
\arguments{
\item{Y}{The matrix containing rows corresponding the curves}

\item{m_not_vector}{The vector containing m_not values for each cluster}

\item{phi_matrix}{A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B}

\item{K}{The total number of clusters}

\item{A_vector}{A vector}

\item{R_vector}{A vector}

\item{probability_matrix}{A matrix in which the rows represent the probabilities that the curve is in each of the clusters}

\item{x}{The x used to generate the clusters}

\item{nbasis}{The number of basis functions}

\item{sigma_list}{A list of the updated sigma parameters for each cluster}
}
\value{
A list of the updated m parameters for each cluster
}
\description{
Update the m parameter for each cluster
}
\examples{
update_m_list(K, A_vector, R_vector, probability_matrix, x, nbasis)
}
