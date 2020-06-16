#' Fit the model to the data
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param function_data The data generated from the actual functions
#' @param K The number of clusters in the data
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @export
#'
#' @return The fitted model
#'
#' @examples fit(Y, K)

fit <- function(Y, K, nbasis, x, init, true_cluster_assignments, gamma_dist_config_matrix, convergence_threshold) {
  print(Y)
}

