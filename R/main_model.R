library(fda)
library(MASS)
source("R/elbo_convergence.R")
source("R/plot.R")

#' Generates cluster assignments and related information given functional data. 
#'
#' @param Y A matrix in which the rows represent the curves 
#' @param K The number of clusters in the data
#' @param nbasis The number of basis functions
#' @param x The dependent variable 
#' @param init The initialization method for the algorithim
#' @param true_cluster_assignments The true cluster assignments 
#' @param gamma_dist_config_matrix A matrix where the rows are the a_not and r_not for each cluster
#' @param d_not_vector A vector of d which is the paramter of Dirichlet distribution for pi
#' @param m_not_vector A matrix of mean vector of phi, the coefficient vector,for all clusters
#' @param v_not_vector A vector of the precision of phi
#' @param convergence_threshold The threshold that determines when the model has converged 
#' @param max_iterations The maximum amount of iterations for the algorithim
#' @param verbose A boolean indicating whether or not to print the inner parameters of the model
#' @param plot_params List of parameters corresponding to characteristics of the plot. Must include vectors xlim and ylim corresponding to the x and y limits of the plot. 
#' @param draw Parameter that determines whether or not to plot the real vs estimated curves after fitting the model 
#' 
#' @return A list with entries containing varius information about the fitted model
#' 
#' @export
#'
#' @examples fit(x, Y, K, true_cluster_assignments, init, nbasis, convergence_threshold, gamma_dist_config_matrix, verbose, draw, plot_params)

funcslustVI <- function(x, Y, K, true_cluster_assignments, init, nbasis, 
                        convergence_threshold, max_iterations, 
                        gamma_dist_config_matrix = NULL, d_not_vector = NULL, m_not_vector = NULL, v_not_vector = NULL,
                        verbose, draw, plot_params) {
  probability_matrix = NULL
  
  if (init == 'hcl') {
    probability_matrix = get_approx_probability_matrix_hcl(Y, K)
  } else if (init == 'tpm') {
    probability_matrix = get_true_probability_matrix(true_cluster_assignments, K)
  } else if (init == "cust") {
    probability_matrix = init
  } else {
    probability_matrix = get_approx_probability_matrix_km(Y, K, x)
  }
  
  init_cluster_assignments = get_final_cluster_assignments(probability_matrix) 
  
  tau_list = get_tau_list(Y, probability_matrix, K) # this can be moved to line 63, else....
  
  if (is.null(true_cluster_assignments) == FALSE) {
    true_m_not = get_true_m_not(x, Y, K, nbasis, true_cluster_assignments)
  }
  if(is.null(m_not_vector)){
    phi_matrix = get_approx_phi_matrix(Y, K, nbasis, probability_matrix, x)
    m_not_vector = phi_matrix
    warning("No prior information for m_not, approximated values have been used based on your data!")
  }
  A_vector = NULL
  R_vector = NULL
  alpha_vector = NULL
  beta_vector = NULL
  if (is.null(gamma_dist_config_matrix) != TRUE) {
    alpha_vector = gamma_dist_config_matrix[1, ]
    observations_per_curve = NCOL(Y)
    A_vector = get_A_vector(alpha_vector, observations_per_curve)
    R_vector = c(gamma_dist_config_matrix[2, ])
    beta_vector = c(gamma_dist_config_matrix[2, ])
  } else {
    alpha_vector = get_alpha_vector(tau_list)
    beta_vector = get_beta_vector(tau_list)
    observations_per_curve = NCOL(Y)
    A_vector = get_A_vector(alpha_vector, observations_per_curve)
    # R_vector = get_beta_vector(tau_list)
    R_vector = beta_vector
    warning("No prior information for a_not and r_not, approximated values have been used based on your data!")
  }
  
  if (verbose == TRUE) {
    print("Probability Matrix Initialization")
    print(probability_matrix)
  }
  
  B = get_B(x, nbasis)
  
  converged = FALSE
  prev_elbo = NULL
  iteration = 0
  while (converged == FALSE & iteration <= max_iterations) {
    iteration = iteration + 1 
    sigma_list = update_sigma_list(Y, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, v_not_vector)
    m_list = update_m_list(Y, m_not_vector, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, sigma_list, v_not_vector)
    R_vector = update_R_vector(Y, K, probability_matrix, x, sigma_list, m_list, nbasis, beta_vector)
    d_vector = update_d_vector(Y, K, probability_matrix, d_not_vector)
    probability_matrix = update_probability_matrix(Y, K, x, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis)
    phi_matrix = get_approx_phi_matrix(Y, K, nbasis, probability_matrix, x)
    elbo = get_elbo(x, Y, K, phi_matrix, m_not_vector, nbasis, sigma_list, m_list, A_vector, R_vector, d_vector, probability_matrix, alpha_vector, beta_vector)
    curr_elbo = elbo[1] 
    converged = check_convergence(prev_elbo, curr_elbo, convergence_threshold)
    prev_elbo = curr_elbo
    
    if (verbose == TRUE) {
      cat("Iteration: ", toString(iteration))
      print("Sigma List")
      print(sigma_list)
      print("M List")
      print(m_list)
      print("R Vector")
      print(R_vector)
      print("D vector")
      print(d_vector)
      print("Probability Matrix")
      print(probability_matrix)
      cat("ELBO: ", toString(curr_elbo))
    }
  }
  # if(itert)
  cluster_assignments = get_final_cluster_assignments(probability_matrix)
  
  if (draw == TRUE & is.null(true_cluster_assignments) == FALSE) {
    plot_data(x, Y, B, m_list, true_m_not, true_cluster_assignments, plot_params)
  }
  result_list = list("probability_matrix" = probability_matrix, "cluster_assignments" = cluster_assignments, 
                     "init_cluster_assignments" = init_cluster_assignments,
                     "Total iteration" = iteration,
                     "A vector" = A_vector,
                     "R vector" = R_vector,
                     "Sigma List" = sigma_list,
                     "M list" = m_list,
                     "Converged ELBO" = elbo,
                     "True basis coefficients" = true_m_not)
  return(result_list)
}

#' Update the sigma parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param K The total number of clusters
#' @param A_vector A vector
#' @param R_vector A vector
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return A list of the updated sigma parameters for each cluster
#'
#' @examples
#' update_sigma_list(Y, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis) 
#'

update_sigma_list <- function(Y, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, v_not_vector) {
  ev_q_tau = A_vector / R_vector
  I = diag(nbasis)
  if(is.null(v_not_vector)){
    v_not_vector = get_v_not_vector(Y, x, probability_matrix, phi_matrix, K, nbasis)
    warning("No prior information for v_not_vector, approximated values have been used based on your data!")
  }
  B = get_B(x, nbasis)
  sigma_list = list()
  for (i in 1:K) {
    sum_matrix = matrix(0, nbasis, nbasis)
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      # v_not = v_not_vector[i]
      temp = p*(t(B) %*% B) ### check, change
      sum_matrix = sum_matrix + temp
    }
    
    mat = ev_q_tau[i] * sum_matrix + v_not_vector[i] * I # here change
    sigma = MASS::ginv(mat)
    sigma_list[[i]] = sigma
  }
  return(sigma_list)
}

#' Update the m parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param m_not_vector The vector containing m_not values for each cluster
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param K The total number of clusters
#' @param A_vector A vector
#' @param R_vector A vector
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#' @param sigma_list A list of the updated sigma parameters for each cluster
#'
#' @return A list of the updated m parameters for each cluster
#'
#' @examples update_m_list(K, A_vector, R_vector, probability_matrix, x, nbasis)

update_m_list <- function(Y, m_not_vector, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, sigma_list, v_not_vector) {
  ev_q_tau = A_vector / R_vector
  if(is.null(v_not_vector)){
    v_not_vector = get_v_not_vector(Y, x, probability_matrix, phi_matrix, K, nbasis)
    warning("No prior information for v_not_vector, approximated values have been used based on your data!")
  }
  B = get_B(x, nbasis)
  m_list = list()
  
  for (i in 1:K) {
    sum_matrix = matrix(0, 1, nbasis)
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      # v_not = v_not_vector[i]
      temp = p * (Y[j, ] %*% B)
      # temp = p * (Y[j, ] %*% B + v_not * m_not_vector[i, ]) # problem
      sum_matrix = sum_matrix + temp
    }
    
    m = (ev_q_tau[i] * sum_matrix + v_not_vector[i] * m_not_vector[i, ]) %*% sigma_list[[i]] # a big problem, fixed
    m_list[[i]] = m
  }
  return(m_list)
}

#' Update the d parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param sigma_list A list of the updated sigma_list parameters for each cluster
#' @param m_list A vlist of the updated m parameters for each cluster
#' @param nbasis The number of basis functions]
#' @param beta_vector A vector of the beta parameters for each cluster
#' 
#' @return A vector of the updated R parameters for each cluster
#'
#' @examples
#' update_R_vector(Y, K, probability_matrix, x, sigma_list, m_list, nbasis, beta_vector)
#'
update_R_vector <- function(Y, K, probability_matrix, x, sigma_list, m_list, nbasis, beta_vector) {
  # tau_list = get_tau_list(Y, probability_matrix, K)
  # beta_vector = get_beta_vector(tau_list, observations_per_curve)
  
  B = get_B(x, nbasis)
  
  R_vector = c(1:K)
  for (i in 1:K) {
    r_not = beta_vector[i]
    sum = 0
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      ev_phi = sum(diag(B %*% sigma_list[[i]] %*% t(B))) + t((Y[j, ] - B %*% t(m_list[[i]]))) %*% (Y[j, ] - B %*% t(m_list[[i]]))
      sum = sum + (p * ev_phi)
    }
    R = r_not + 1/2 * sum
    R_vector[i] = R
  }
  return(R_vector)
}

#' Update the d parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' 
#' @return A vector of the updated d parameters for each cluster
#'
#' @examples
#' update_d_vector(Y, K, probability_matrix)

update_d_vector <- function(Y, K, probability_matrix, d_not_vector) {
  if(is.null(d_not_vector)){
    d_not_vector = get_d_not_vector(K)
    warning("No prior information for d_not, 1/K has been used for each cluster.")
  }
  
  d_vector = c(1:K)*0
  for (i in 1:K) {
    sum = 0
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      sum = sum + p
    }
    d = d_not_vector[i] + sum
    d_vector[i] = d
  }
  return(d_vector)
}

# update_d_vector <- function(K, probability_matrix) {
# d_not_vector = get_d_not_vector(K)
# d = apply(probability_matrix, 2, sum) + d_not_vector
# return(d)
# }

#' Update the probabilty matrix
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param m_list A list of the m parameters for each cluster
#' @param A_vector A vector
#' @param R_vector A vector
#' @param d_vector A vector of the d parameters for each cluster
#' @param nbasis The number of basis functions
#' 
#' @return The updated probability matrix 
#' 
#' @examples
#' update_probability_matrix(Y, K, probability_matrix, x, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis)

update_probability_matrix <- function(Y, K, x, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis) {
  probability_matrix = matrix(0, NROW(Y), K)
  observations_per_curve = NCOL(Y)
  B = get_B(x, nbasis)
  
  coef = observations_per_curve / 2
  for (k in 1:K) {
    for (i in 1:NROW(Y)) {
      ev_phi = sum(diag(B %*% sigma_list[[k]] %*% t(B))) + t((Y[i, ] - B %*% t(m_list[[k]]))) %*% (Y[i, ] - B %*% t(m_list[[k]]))
      ev_tau_k_log_tau_k = digamma(A_vector[k]) - log(R_vector[k]) # same problem fixed, not log10 but log
      ev_tau_k_tau_k = A_vector[k] / R_vector[k]
      ev_pi_log_pi_k = digamma(d_vector[k]) - digamma(sum((d_vector)))
      
      a_i_k = coef * ev_tau_k_log_tau_k - 1/2 * ev_tau_k_tau_k * ev_phi + ev_pi_log_pi_k
      num = exp(a_i_k)
      
      # sum = 0
      # for (k2 in 1:K) {
      # ev_phi = sum(diag(B %*% sigma_list[[k2]] %*% t(B))) + t((Y[i, ] - B %*% t(m_list[[k2]]))) %*% (Y[i, ] - B %*% t(m_list[[k2]]))
      #  ev_tau_k_log_tau_k = digamma(A_vector[k2]) - log10(R_vector[k2])
      # ev_tau_k_tau_k = A_vector[k2] / R_vector[k2]
      # ev_pi_log_tau_k = digamma(d_vector[k2]) - digamma(sum((d_vector)))
      # sum_a_i_k = coef * ev_tau_k_log_tau_k - 1/2 * ev_tau_k_tau_k * ev_phi + ev_pi_log_tau_k
      #  sum = sum + exp(sum_a_i_k)
      # }
      
      # if ((num == 0) & (sum == 0)) {
      # p = 0
      # } else {
      # p = round(num / sum, 4)
      # }
      probability_matrix[i, k] = num
      # probability_matrix[i, k] = p
    }
  }
  probability_matrix_row_sum = apply(probability_matrix, 1, sum)
  for (i in 1:nrow(probability_matrix)) {
    if(probability_matrix_row_sum[i] == 0) {
      maxindex <- which.max(probability_matrix[i, ])
      probability_matrix[i, maxindex] <- 1
      probability_matrix[i, -maxindex] <- 0
    }
    # probability_matrix[i, ] <- 0
    else probability_matrix[i, ] <- probability_matrix[i, ] / probability_matrix_row_sum[i]
  }
  return(probability_matrix)
}

#' Generates a matrix in which each row represents the variances of the curves in each cluster
#'
#' @param Y The matrix containing rows corresponding to the curves
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param K The number of entries in each curve vector
#'
#' @return A matrix in which each row represents the variances of the curves in each cluster
#'
#' @examples
#' get_tau_list(Y, probability_matrix, K)

get_tau_list <- function(Y, probability_matrix, K) {
  cumulative_matrix_list = vector("list", length = K)
  
  for (i in 1:NROW(Y)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    if(is.null(cumulative_matrix_list[[max_prob_index]]) == TRUE) {
      mat = matrix(0, 1, NCOL(Y))
      mat[1, ] = Y[i, ]
      cumulative_matrix_list[[max_prob_index]] = mat
    } else {
      cumulative_matrix_list[[max_prob_index]] = rbind(cumulative_matrix_list[[max_prob_index]], Y[i, ])
    }
  }
  # the purpose of this part above is to group the curves into K clusters
  tau_list = vector("list", length = K)
  # Var_list = vector("list", length = K)
  for (i in 1:K) {
    vec <- vector()
    tau_list[[i]] <- vec
    # Var_list[[i]] <- vec
  }
  
  for (i in 1:K) {
    mat = cumulative_matrix_list[[i]]
    
    if (is.null(mat) == TRUE) {
      tau_list[[i]] = c(1:NCOL(Y))*0 # I think here should make some changes, set to a vey small value?
    } else {
      for(j in 1:NCOL(Y)) {
        Var = var(mat[,j])
        # change made here, Sept 19
        if(is.na(Var) == TRUE) Var = 10^(-20)
        tau = 1 / Var # change var  to Var to avoid repeat use of default function in R
        tau_list[[i]] = c(tau_list[[i]], tau) 
        # Var_list[[i]] = c(Var_list[[i]], Var)
      }
    }
    # tau_list[[i]] <- rep(1 / mean(Var_list[[i]]), NCOL(Y)) # change made here, Sept 27 
  }
  
  return(tau_list)
  
}

#' Generates a vector A
#'
#' @param alpha_vector A vector that in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#' @param observations_per_curve The number of entries in each curve vector
#'
#' @return A vector A
#'
#' @examples
#' get_A(alpha_vector, observations_per_curv)

get_A_vector <- function(alpha_vector, observations_per_curve) {
  A_vector = alpha_vector + observations_per_curve / 2
  return(A_vector)
}

#' Generates a vector in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @param tau_list A list in which each index is the tau vector for that cluster
#'
#' @return A vector that in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @examples
#' get_alpha_vector(tau_list)

get_alpha_vector <- function(tau_list) {
  alpha_vector = c(1:length(tau_list))*0
  for (cluster_number in 1:length(tau_list)) {
    tau_k = tau_list[[cluster_number]]
    # should consider variance = 0, that's  tau_k are the same through different time points in cluster k
    if(length(unique(tau_k)) == 1) {
      alpha_vector[cluster_number] = 0 
    } else {
      expected_value = mean(tau_k)
      variance = var(tau_k)
      alpha = expected_value ^ 2 / variance # use Gamma property
      alpha_vector[cluster_number] = alpha
    }
  }
  return(alpha_vector)
}

#' Generates a vector in which the entries are the beta parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @param tau_list A list in which each index is the tau vector for that cluster
#'
#' @return A vector that in which the entries are the beta parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @examples
#' get_beta_vector(tau_list, observations_per_cluster_x)

get_beta_vector <- function(tau_list) {
  beta_vector = c(1:length(tau_list))*0
  for (cluster_number in 1:length(tau_list)) {
    tau_k = tau_list[[cluster_number]]
    if(length(unique(tau_k)) == 1) {
      beta_vector[cluster_number] = 0 # set as 0?
    } else {
      expected_value = mean(tau_k)
      variance = var(tau_k)
      beta = expected_value / variance
      beta_vector[cluster_number] = beta
    }
  }
  return(beta_vector)
}

#' Generates initial probability matrix with approximated probability values via K-mean clustering for each cluster
#' 
#' @param Y The matrix containing rows corresponding the curves
#' @param K The number of clusters in cluster data
#' @param x The dependent variable 
#'
#' @return A probability matrix with approximated probability values for each cluster using kmeans 
#'
#' @examples get_approx_probability_matrix(Y, K, x)

get_approx_probability_matrix_km <- function(Y, K, x) {  # In this function, x is not necessary currently
  res = kmeans(Y, K) 
  predictions = res$cluster
  probability_matrix = matrix(0, NROW(Y), K)
  for (i in 1:length(predictions)) {
    cluster_prediction = predictions[[i]]
    probability_matrix[i, cluster_prediction] = 1
  }
  
  return(probability_matrix)
}

#' Generates initial probability matrix with approximated probability values via hierarchical clustering for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The number of clusters in cluster data
#' @param x The x used to generate the clusters
#'
#' @return A probability matrix with probability values for each cluster
#'
#' @examples get_equal_probability_matrix(K, curves_per_cluster, x)

get_approx_probability_matrix_hcl <- function(Y, K, x) {
  predictions = cutree(hclust(dist(Y), method = "ward"), K) # use distance for clustering
  probability_matrix = matrix(0, NROW(Y), K)
  for (i in 1:length(predictions)) {
    cluster_prediction = predictions[[i]]
    probability_matrix[i, cluster_prediction] = 1
  }
  
  return(probability_matrix)
}

#' Generates B matrix
#'
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return A matrix in which row represent curves of various clusters
#'
#' @examples
#' get_B(x, nbasis)

get_B <- function(x, nbasis) {
  rangeval = c(0, x[length(x)])
  basisobj = fda::create.bspline.basis(rangeval, nbasis)
  B <- fda::getbasismatrix(x, basisobj=basisobj)
  return(B)
}

#' Gets matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param nbasis The number of basis functions
#' @param probability_matrix The x used to generate the clusters
#' @param x The x used to generate the clusters
#' 
#' @return A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#'
#' @examples
#' get_approx_phi_matrix(Y, K, nbasis, probability_matrix, x)

get_approx_phi_matrix <- function(Y, K, nbasis, probability_matrix, x) {
  function_data = get_approx_function_data(Y, K, probability_matrix)
  observations_per_curve = length(x)
  min_arg = x[1]
  max_arg = x[observations_per_curve]
  basisobj = fda::create.bspline.basis(c(min_arg, max_arg), nbasis)
  phi_matrix = matrix(0, K, nbasis)
  for (i in 1:K) {
    f = function_data[i, ]
    phi = fda::smooth.basis(argvals = x, y = f, fdParobj = basisobj)$fd$coef
    phi_matrix[i, ] = c(phi)
  }
  return(phi_matrix)
}

#' Gets the vnot parameter
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param x The x used to generate the clusters
#' @param probability_matrix The x used to generate the clusters
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param K The total number of clusters
#' @param nbasis The number of basis functions
#'
#' @return The vnot parameter
#'
#' @examples
#' get_vnot(Y, x, probability_matrix, phi_matrix, K, nbasis)

get_v_not_vector <- function(Y, x, probability_matrix, phi_matrix, K, nbasis) {
  cumulative_phi_matrix = get_cumulative_phi_matrix(Y, x, nbasis)
  cumulative_phi_matrix_list = vector("list", length = K)
  
  for (i in 1:NROW(Y)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    if (is.null(cumulative_phi_matrix_list[[max_prob_index]]) == TRUE)   {
      mat = matrix(0, 1, nbasis)
      mat[1, ] = cumulative_phi_matrix[i, ]
      cumulative_phi_matrix_list[[max_prob_index]] = mat
    } else {
      cumulative_phi_matrix_list[[max_prob_index]] = rbind(cumulative_phi_matrix_list[[max_prob_index]], cumulative_phi_matrix[i, ])
    }
  }
  
  v_not_vector = c(1:K)
  for(i in 1:K) {
    var_sum = 0
    mat = cumulative_phi_matrix_list[[i]]
    if (is.null(mat) == TRUE) {
      v_not = 0
    } else {
      for(j in 1:nbasis) {
        var_sum = var_sum + var(mat[, j])
      } 
      # change made here Sept 19
      if(is.na(var_sum) == TRUE) var_sum = 10^(-20)
      avg_var = var_sum / nbasis
      v_not = 1 / avg_var
    }
    v_not_vector[i] = v_not
  }
  return(v_not_vector)
}

#' Gets the d_not parameter for each cluster
#'
#' @param K The total number of clusters
#'
#' @return The d_not parameter vector for each cluster
#'
#' @examples
#' get_d_not_vector(K)

get_d_not_vector <- function(K) {
  val = 1 / K
  d_not_vector = c(1:K)
  
  for (i in 1:K) {
    d_not_vector[i] = val
  }
  
  return(d_not_vector)
}

#' Gets functional data approximations, use mean for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#'
#' @return Matrix with approximations of functional data for each cluster in a row
#'
#' @examples
#' get_approx_function_data(Y, number_of_cluster)

get_approx_function_data <- function(Y, K, probability_matrix) {
  sum_count_vector = c(1:K)*0
  sum_matrix = matrix(0, K, NCOL(Y))
  
  for (i in 1:NROW(Y)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    sum_count_vector[max_prob_index] = sum_count_vector[max_prob_index] + 1
    sum_matrix[max_prob_index, ] = sum_matrix[max_prob_index, ] + Y[i, ]
  }
  
  function_data = matrix(0, K, NCOL(Y))
  for (i in 1:K) {
    function_data[i, ] = sum_matrix[i, ] / sum_count_vector[i]
  }
  
  return(function_data)
}

#' Gets a matrix of the coffecient vectors, phi, for each curve
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return a matrix of the coffecient vectors, phi, for each curve
#'
#' @examples
#' get_cumulative_phi_matrix(Y, x, nbasis)

get_cumulative_phi_matrix <- function(Y, x, nbasis) {
  observations_per_curve = length(x)
  min_arg = x[1]
  max_arg = x[observations_per_curve]
  basisobj = fda::create.bspline.basis(c(min_arg, max_arg), nbasis)
  cumulative_phi_matrix = matrix(0, NROW(Y), nbasis)
  for (i in 1:NROW(Y)) {
    f = Y[i, ]
    phi = fda::smooth.basis(argvals = x, y = f, fdParobj = basisobj)$fd$coef
    cumulative_phi_matrix[i, ] = c(phi)
  }
  return(cumulative_phi_matrix)
}

#' Gets a vector of the final cluster assignments based on the probability matrix 
#'
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' 
#' @return A vector of the final cluster assignments
#'
#' @examples
#' get_final_cluster_assignments(probability_matrix)

get_final_cluster_assignments <- function(probability_matrix) {
  final_cluster_assignments = c(1:NROW(probability_matrix))
  
  for (i in 1:NROW(probability_matrix)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    final_cluster_assignments[i] = max_prob_index
  }
  return(final_cluster_assignments)
}

#' Gets the true probability matrix using the true cluster assignments 
#'
#' @param true_cluster_assignments A vector containing the true cluster assignments 
#' @param K The number of clusters in the data 
#' 
#' @return The true probability matrix of the true cluster assignments
#'
#' @examples
#' get_true_probability_matrix(true_cluster_assignments, K)

get_true_probability_matrix <- function(true_cluster_assignments, K) {
  number_of_curves = length(true_cluster_assignments)
  probability_matrix = matrix(0, number_of_curves, K)
  for (i in 1:number_of_curves) {
    col = true_cluster_assignments[i]
    probability_matrix[i, col] = 1
  }
  return(probability_matrix)
}

#' Gets the vector with the true m_not values
#'
#' @param x The x used to generate the clusters
#' @param Y The matrix containing rows corresponding the curves
#' @param K The number of clusters in the data 
#' @param nbasis The number of basis functions
#' @param true_cluster_assignments A vector containing the true cluster assignments 
#' 
#' 
#' @return The vector with the true m_not values 
#'
#' @examples
#' get_true_m_no(x, Y, K, nbasis, true_cluster_assignments)


get_true_m_not <- function(x, Y, K, nbasis, true_cluster_assignments) {
  true_probability_matrix = get_true_probability_matrix(true_cluster_assignments, K)
  true_m_not = get_approx_phi_matrix(Y, K, nbasis, true_probability_matrix, x)
  return(true_m_not)
}

#' @param m_list A list of the updated m parameters for each cluster
#' @param true_m_not A matrix containing the true m_not vectors for each cluster

get_mse <- function(B, true_m_not, m_list){
  number_of_clusters = NROW(true_m_not)
  mse <- numeric(number_of_clusters)
  for (i in 1:number_of_clusters) {
    mse[i] <- mean(((B %*% true_m_not[i, ]) - (B %*% t(m_list[[i]])))^2) 
  }
  return(mse)
}


