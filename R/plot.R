#' Generates a plot with the true vs predicted curves
#'
#' @param x The dependent variable 
#' @param Y A matrix in which the rows represent the curves 
#' @param K The number of clusters in the data
#' @param nbasis The number of basis functions
#' @param m_list A list of the updated m parameters for each cluster
#' @param true_m_not A matrix containing the true m_not vectors for each cluster
#'
#' @examples plot_data(x, Y, K, nbasis, m_list, true_m_not)

plot_data <- function(x, B, m_list, true_m_not) {
  number_of_clusters = NROW(true_m_not)
  
  plot(x, B %*% true_m_not[1, ], col=1, lwd=2, type="l", ylim=c(1, 6), main="Plot", ylab="f(x)", xlab="x")
  lines(x, B %*% t(m_list[[1]]), col=2, lwd=2)
  
  for (i in 2:number_of_clusters) {
    lines(x, B %*% true_m_not[i, ], col=1, lwd=2)
    lines(x, B %*% t(m_list[[i]]), col=2, lwd=2)
  }
}
