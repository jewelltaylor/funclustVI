
#' Gets simulated data from Case_7 
#'
#' @param x The x used to generate the clusters
#' @param n The number of curves per cluster 
#' 
#' @return The simulated data from Case_7 
#' 
#' @export
#'
#' @examples
#' Case_7(x, n)

Case_7 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/4,1/4) + .3+ (1/1.3)*sin(x*1.3) + x^3 + rnorm(length(x),2,0.4^2)
    Y[,i+n] = runif(1,-1/4,1/4) + 1.0+ (1/1.2)*sin(x*1.3) + x^3 + rnorm(length(x),2,0.4^2)
    Y[,i+2*n] = runif(1,-1/4,1/4) + .2+ (1/4)*sin(x*1.3) + x^3 + rnorm(length(x),2,0.4^2)
  }
  return(t(Y))
}

#' Gets the simulated data from Case_8 
#'
#' @param x The x used to generate the clusters
#' @param n The number of curves per cluster 
#' 
#' @return The simulated data from Case_8
#' 
#' @export
#'
#' @examples
#' Case_8(x, n)

Case_8 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/2,1/2) +1.1+ sin(pi*x*1.5) + cos(pi*x^2) + rnorm(length(x),2,0.4^2)
    Y[,i+n] = runif(1,-1/2,1/2) +1.5+ sin(pi*x*1.7) + cos(pi*x^2) + rnorm(length(x),2,0.4^2)
    Y[,i+2*n] = runif(1,-1/2,1/2) +2.2+ sin(pi*x*1.9) + cos(pi*x^2) + rnorm(length(x),2,0.4^2)
  }
  return(t(Y))
}

#' Gets the simulated data from Case_9
#'
#' @param x The x used to generate the clusters
#' @param n The number of curves per cluster 
#' 
#' @return The simulated data from Case_9
#' 
#' @export
#'
#' @examples
#' Case_9(x, n)

Case_9 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/4,1/4) + (1/1.8)*exp(x*1.1) - x^3+ rnorm(length(x),2,0.3^2)
    Y[,i+n] = runif(1,-1/4,1/4) + (1/1.7)*exp(x*1.4) - x^3+ rnorm(length(x),2,0.3^2)
    Y[,i+2*n] = runif(1,-1/4,1/4) + (1/1.5)*exp(x*1.5) - x^3+ rnorm(length(x),2,0.3^2)
  }
  return(t(Y))
}

#' Gets the simulated data from Case_11
#'
#' @param x The x used to generate the clusters
#' @param n The number of curves per cluster 
#' 
#' @return The simulated data from Case_11
#' 
#' @export
#'
#' @examples
#' Case_11(x, n)


Case_11 <- function(data_params) # 4 clusters
{
  K = 4
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/3,1/3) +0.2- sin(pi*x*1.1) +x^3 + rnorm(length(x),2,0.4^2)
    Y[,i+n] = runif(1,-1/3,1/3) +0.5- sin(pi*x*1.4) +x^3 + rnorm(length(x),2,0.4^2)
    Y[,i+2*n] = runif(1,-1/3,1/3) +0.7- sin(pi*x*1.6) +x^3 + rnorm(length(x),2,0.4^2)
    Y[,i+3*n] = runif(1,-1/3,1/3) +1.3- sin(pi*x*1.8) +x^3 + rnorm(length(x),2,0.4^2)
  }
  return(t(Y))
}

#' Gets the simulated data from Case_12
#'
#' @param x The x used to generate the clusters
#' @param n The number of curves per cluster 
#' 
#' @return The simulated data from Case_12
#' 
#' @export
#'
#' @examples
#' Case_12(x, n)

Case_12 <- function(data_params) # 5 clusters
{
  K = 5
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/3,1/3) +1.1+ sin(pi*x*1.5) + sin(pi*x^2) + rnorm(length(x),2,0.3^2)
    Y[,i+n] = runif(1,-1/3,1/3) +1.5+ sin(pi*x*1.7) + sin(pi*x^2) + rnorm(length(x),2,0.3^2)
    Y[,i+2*n] = runif(1,-1/3,1/3) +2.2+ sin(pi*x*1.9) + sin(pi*x^2) + rnorm(length(x),2,0.3^2)
    Y[,i+3*n] = runif(1,-1/3,1/3) +2.3+ sin(pi*x*1.8) + sin(pi*x^2) + rnorm(length(x),2,0.3^2)
    Y[,i+4*n] = runif(1,-1/3,1/3) +2.4+ sin(pi*x*1.6) + sin(pi*x^2) + rnorm(length(x),2,0.3^2)
  }
  return(Y)
}

#' Gets the simulated data from Case_44
#'
#' @param x The x used to generate the clusters
#' @param n The number of curves per cluster 
#' 
#' @return The simulated data from Case_44
#' 
#' @export
#'
#' @examples
#' Case_44(x, n)

Case_44 <- function(data_params) # 6 clusters
{
  K = 6
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/4,1/4) + cos(pi*x*1) - x^2 + rnorm(length(x),2,0.3^2)
    Y[,i+n] = runif(1,-1/4,1/4) + cos(pi*x*1.2) - x^2 + rnorm(length(x),2,0.3^2)
    Y[,i+2*n] = runif(1,-1/4,1/4) + cos(pi*x*1.4) - x^2 + rnorm(length(x),2,0.3^2)
    Y[,i+3*n] = runif(1,-1/4,1/4) + cos(pi*x*1.6) - x^2 + rnorm(length(x),2,0.3^2)
    Y[,i+4*n] = runif(1,-1/4,1/4) + cos(pi*x*1.8) - x^2 + rnorm(length(x),2,0.3^2)
    Y[,i+5*n] = runif(1,-1/4,1/4) + cos(pi*x*2) - x^2 + rnorm(length(x),2,0.3^2)
  }
  return(t(Y))
}
