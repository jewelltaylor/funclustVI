library(fda)

#' Gets the elbo value for a single iteration of the algorithim 
#'
#' @param x The x used to generate the clusters
#' @param Y The matrix containing rows corresponding the curves
#' @param K The number of clusters in the data 
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param m_not_vector The vector containing m_not values for each cluster
#' @param nbasis The number of basis functions
#' @param m_list A list of the m parameters for each cluster
#' @param A_vector A vector
#' @param R_vector A vector
#' @param d_vector A vector of the d parameters for each cluster
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param alpha_vector A vector that in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#' @param beta_vector A vector contain the beta vector for each cluster
#' 
#' @return The elbo value 
#'
#' @examples
#' get_elbo(x, Y, K, phi_matrix, m_not_vector, nbasis, sigma_list, m_list, A_vector, R_vector, d_vector, probability_matrix, alpha_vector, beta_vector)

get_elbo <- function(x, Y, K, phi_matrix, m_not_vector, nbasis, sigma_list, m_list, A_vector, R_vector, d_vector, probability_matrix, alpha_vector, beta_vector) {
  B = get_B(x, nbasis)
  d_not_vector = get_d_not_vector(K)
  diff_pi <- sum((d_not_vector - d_vector) * (digamma(d_vector) - digamma(sum(d_vector)))) 
  
  diff_z_sum_1 = 0
  diff_z_sum_2 = 0
  for (i in 1:NROW(Y)) {
    for (k in 1:K) {
      ev_pi_log_tau_k = digamma(d_vector[k]) - digamma(sum((d_vector)))
      diff_z_sum_1 = probability_matrix[i, k] * ev_pi_log_tau_k + diff_z_sum_1
      if (probability_matrix[i, k] <= .0000001) {
        diff_z_sum_2 = 0 + diff_z_sum_2
      } else {
        diff_z_sum_2 = probability_matrix[i, k] * log(probability_matrix[i, k]) + diff_z_sum_2
      }
    }
  }
  diff_z = diff_z_sum_1 - diff_z_sum_2
  
  v_not_vector <- get_v_not_vector(Y, x, probability_matrix, phi_matrix, K, nbasis)
  diff_phi_sum <- 0
  for (k in 1:K) {
    diff_phi_sum <- - sum((v_not_vector[k] * (sum(diag(sigma_list[[k]])) + (m_list[[k]] - m_not_vector[k, ]) %*% 
                                                t(m_list[[k]] - m_not_vector[k, ])) + log(det(sigma_list[[k]])))) / 2 + diff_phi_sum
  }
  diff_phi <- diff_phi_sum
  
  tau_list = get_tau_list(Y, probability_matrix, K)
  diff_tau_sum_1 = 0
  diff_tau_sum_2 = 0
  for (k in 1:K) {
    ev_tau_k_log_tau_k = digamma(A_vector[k]) - log(R_vector[k])
    diff_tau_sum_1 = (alpha_vector[k] - 1) * ev_tau_k_log_tau_k - beta_vector[k] * (A_vector[k] / R_vector[k]) + diff_tau_sum_1
    diff_tau_sum_2 = A_vector[k] * log(R_vector[k]) + (A_vector[k] - 1) * ev_tau_k_log_tau_k + diff_tau_sum_2
  }
  diff_tau = diff_tau_sum_1 - diff_tau_sum_2
  
  ev_log_lik_sum = 0
  for (i in 1:NROW(Y)) {
    for (k in 1:K) {
      ev_tau_k_log_tau_k = digamma(A_vector[k]) - log(R_vector[k])
      ev_phi = sum(diag(B %*% sigma_list[[k]] %*% t(B))) + t((Y[i, ] - B %*% t(m_list[[k]]))) %*% (Y[i, ] - B %*% t(m_list[[k]]))
      ev_log_lik_sum = probability_matrix[i, k] * (1/2 * ev_tau_k_log_tau_k - 1/2 * (A_vector[k] / R_vector[k]) * ev_phi) + ev_log_lik_sum
    }
  }
  ev_log_lik = as.numeric(ev_log_lik_sum)
  
  elbo = diff_pi + diff_z + diff_phi + diff_tau + ev_log_lik
  return(c(elbo, ev_log_lik))
}
#' Checks if the algorithim has converged with the given threshold 
#'
#' @param prev_elbo The elbo value from the previous iteration 
#' @param curr_elbo The elbo value from the current iteration 
#' @param convergence_threshold The threshold that determines when the model has converged 
#' 
#' 
#' @return A boolean whether the algorithim has converged 
#'
#' @examples
#' check_convergence(prev_elbo, curr_elbo, convergence_threshold)

check_convergence <- function(prev_elbo, curr_elbo, convergence_threshold) {
  if(is.null(prev_elbo) == TRUE) {
    return(FALSE)
  }
  else{
    dif = abs(curr_elbo - prev_elbo)
    if(dif  <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
}
