library(gtools)

#' Gets a matrix of the final cluster assignments based on the probability matrix 
#'
#' @param predicited_clusters The predicted clusters 
#' @param true_clusters THe true clusters
#' @param K The number of clusters 
#' 
#' @return a matrix of the final cluster assignments
#'
#' @examples
#' Mismatch(predicited_clusters, true_clusters, K)

Mismatch <- function(predicited_clusters, true_clusters, K)
{
  sigma = gtools::permutations(n=K,r=K,v=1:K)
  
  Miss = length(which( true_clusters != predicited_clusters))  ## for permutation 1, 2,... K
  
  mm_aux = predicited_clusters
  for (ind in 2:dim(sigma)[1])
  {
    
    for (j in 1:K)
      mm_aux[which(predicited_clusters == j)] = sigma[ind,j]
    
    Miss[ind] =  length(which( true_clusters != mm_aux))
    
  }
  return(min(Miss))
}
