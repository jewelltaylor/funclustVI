library(fda)
library(sabre)
source("R/simulate_data.R")
source("R/model.R")
source("R/evaluation_metrics.R")


#' Runs and evaluates the model on the specified simulated function data 
#'
#' @param case The case used to simulate the data
#' @param x The dependent variable 
#' @param K The number of clusters in the data
#' @param curves_per_cluster The number of curves per cluster 
#' @param number_of_simulations The number of simulations to evatuate the model over 
#' @param initialization_method The initialization method for the algorithim
#' @param nbasis The number of basis functions
#' @param true_cluster_assignments The true cluster assignments 
#' @param gamma_dist_config_matrix A matrix where the rows are the alpha and parameters for each cluster
#' @param convergence_threshold The threshold that determines when the model has converged 
#' @param save_path Pathway to save 
#' @param verbose A boolean indicating whether or not to print the inner parameters of the model
#' @param draw A boolean indicating whether or not to draw the true function vs the 
#'
#' @return A matrix with the results of the algorithim for iteration in the simulation for each metric being evaluated over
#' 
#' @export
#'
#' @examples run_simulations(case, x, K, curves_per_cluster, number_of_simulations, initialization_method, nbasis, true_cluster_assignments, gamma_dist_config_matrix, convergence_threshold, save_path, verbose, draw)


run_simulations <- function(data_params, model_params, eval_func_list, number_of_simulations, save_path) {
  num_of_eval_funcs = length(eval_func_list)
  eval_func_name_list = names(eval_func_list)
  eval_metric_sum_vector = c(1:num_of_eval_funcs)*0 
  final_res_mat = matrix(0, number_of_simulations, num_of_eval_funcs)
  
  count = 0
  for (simulation_number in 1:number_of_simulations) {
    res = compute_function(simulation_number, data_params, model_params, eval_func_list)
    
    for (i in 1:num_of_eval_funcs) {
      prev_sum = eval_metric_sum_vector[i]
      curr_val = res[[i]]
      final_res_mat[simulation_number, i] = curr_val
      eval_metric_sum_vector[i] = prev_sum + curr_val
      eval_func_name = eval_func_name_list[i]
      cat(eval_func_name, " = ", curr_val)
    }
    cat("\n")
    count = count + 1
  }
  
  eval_metric_avg_vector = eval_metric_sum_vector / count
  
  for(i in 1:num_of_eval_funcs) {
    eval_func_name = eval_func_name_list[i]
    avg_val = eval_metric_avg_vector[i]
    cat("Average ", eval_func_name, " = ", avg_val)
  }
  cat("\n")
  
  if(is.null(save_path) == FALSE) {
    write(ffinal_res_mat, save_path)
  }
  
  return(final_res_mat)
}

compute_function <- function(simulation_number, data_params, model_params, eval_func_list) {
  set.seed(simulation_number)
  num_of_eval_funcs = length(eval_func_list)
  eval_func_name_list = names(eval_func_list)
  
  generate_data = data_params$generate_data
  
  Y = generate_data(data_params)
  
  eval_metric_res_vector = c(1:num_of_eval_funcs)*0
  
  model_func = model_params$model_func
  
  cluster_assignments = model_func(Y, data_params, model_params)
  
  result_list = list()
  
  for(i in 1:num_of_eval_funcs) {
    eval_func_name = eval_func_name_list[i]
    eval_func = eval_func_list[[i]]
    res = eval_func(cluster_assignments, data_params)
    result_list[[eval_func_name]] = res
  }
  
  return(result_list)
}

get_funclustVI_cluster_assignments <- function(Y, data_params, model_params) {
  x = data_params$x
  K = data_params$K
  init = model_params$init
  nbasis = model_params$nbasis
  convergence_threshold = model_params$convergence_threshold
  gamma_dist_config_matrix = model_params$gamma_dist_config_matrix
  true_cluster_assignments = data_params$true_cluster_assignments
  verbose = model_params$verbose 
  draw = model_params$draw
  clf = funcslustVI(Y, x, K, init, nbasis, convergence_threshold, gamma_dist_config_matrix, true_cluster_assignments, verbose, draw)
  cluster_assignments = clf$cluster_assignments
  return(cluster_assignments)
}


#' Gets the number of mismatches for a single iteration of the simulation 
#'
#' @param cluster_assignments A vector where each entry is the cluster assignment for the corresponging curve 
#' @param K The number of clusters in the data
#' @param curves_per_cluster The number of curves per cluster 
#'
#' @return The number of mismatches  
#'
#' @examples get_mismatches(cluster_assignments, K, curves_per_cluster)


get_mismatches <- function(cluster_assignments, data_params) {
  true_cluster_assignments = data_params$true_cluster_assignments
  mismatches = Mismatch(cluster_assignments, true_cluster_assignments, K)
  return(mismatches)
}

#' Gets the number of vmeaure for a single iteration of the simulation 
#'
#' @param cluster_assignments A vector where each entry is the cluster assignment for the corresponging curve 
#' @param K The number of clusters in the data
#' @param curves_per_cluster The number of curves per cluster 
#'
#' @return The vmeasure
#'
#' @examples get_v_measure(cluster_assignments, K, curves_per_cluster)

get_v_measure <- function(cluster_assignments, data_params) {
  true_cluster_assignments = data_params$true_cluster_assignments
  v_measure = sabre::vmeasure(cluster_assignments, true_cluster_assignments)$v_measure
  return(v_measure)
}


