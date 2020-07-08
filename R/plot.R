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

plot_data <- function(x, Y, B, m_list, true_m_not, true_cluster_assignments, plot_params) {
  number_of_clusters = NROW(true_m_not)
  ylim = plot_params$ylim
  xlim = plot_params$xlim 
  show_curves = plot_params$show_curves 
  title = plot_params$title
  
  if(is.null(title) == TRUE) {
    title = "Plot"
  }
  
  if(isTRUE(show_curves) == TRUE) {
    cluster_number = true_cluster_assignments[1]
    col = 2 + cluster_number
    plot(x, Y[1, ], col=col, type="l", ylim=ylim, xlim=xlim, main=title, ylab="f(x)", xlab="x")
    for(i in 2:NROW(Y)) {
      cluster_number = true_cluster_assignments[i]
      col = 2 + cluster_number
      lines(x, Y[i, ], type="l", col=col)
    }
    
    lines(x, B %*% true_m_not[1, ], col=1, lwd=2)
    lines(x, B %*% t(m_list[[1]]), col=2, lwd=2)
    
    for (i in 2:number_of_clusters) {
      lines(x, B %*% true_m_not[i, ], col=1, lwd=2)
      lines(x, B %*% t(m_list[[i]]), col=2, lwd=2)
    }
    
  } else {
    plot(x, B %*% true_m_not[1, ], col=1, lwd=2, type="l", ylim=ylim, xlim=xlim, main=title, ylab="f(x)", xlab="x")
    lines(x, B %*% t(m_list[[1]]), col=2, lwd=2)
    
    for (i in 2:number_of_clusters) {
      lines(x, B %*% true_m_not[i, ], col=1, lwd=2)
      lines(x, B %*% t(m_list[[i]]), col=2, lwd=2)
    }
  }
  
  legend(x= "topleft", cex=.75, legend=c("True Function", "Estimated Function", "Actual Curves"), col=c("black", "red", "grey"), lty=c(1,1,1))
}
