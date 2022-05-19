## For scenario 1--------------------------------------------------------------------------------------------------------------
# get MSE for each repetition
get_EMSE_1 <- function(repeated_time){
  MSE_list <- list()
  for (i in 1:repeated_time) {
    number_of_simulations = 1
    save_path = NULL
    
    # Data Parameters
    x = seq(from=0,to=pi/3, length = 100)
    K = 3
    curves_per_cluster = 50 
    true_cluster_assignments = rep(1:K,each = curves_per_cluster)
    seeds = i
    
    # Pack into data parameter list
    data_params = list()
    data_params$x = x
    data_params$K = K
    data_params$curves_per_cluster = curves_per_cluster
    data_params$true_cluster_assignments = true_cluster_assignments
    data_params$seeds = seeds 
    data_params$generate_data = Case_7
    
    # Model Parameters
    init = "km"
    nbasis = 6
    gamma_dist_config_matrix = matrix(0, 2, K)
    gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 25
    gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 25
    d_not_vector = c(1/3, 1/3, 1/3)
    m_not_vector = matrix(c(0.19, 0.31, 0.54, 1.08, 1.69, 2.09,
                            1.14, 1.27, 1.52, 2.08, 2.70, 3.11,
                            0.15, 0.19, 0.27, 0.61, 1.16, 1.55), nrow = 3, byrow = T)
    v_not_vector = c(10, 10, 10, 10, 10, 10)
    convergence_threshold = 0.01
    max_iterations = 100
    verbose = FALSE
    draw = FALSE
    
    # Pack into model parameter list 
    model_params = list()
    model_params$model_func = get_funclustVI_cluster_assignments
    model_params$init = "km"
    model_params$nbasis = 6
    model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
    model_params$d_not_vector = d_not_vector
    model_params$m_not_vector = m_not_vector 
    model_params$v_not_vector = v_not_vector   
    model_params$convergence_threshold = convergence_threshold
    model_params$max_iterations = max_iterations 
    model_params$save_path = save_path
    model_params$verbose = verbose
    model_params$draw = draw
    plot_params = list()
    plot_params$xlim = NULL
    plot_params$ylim = c(1, 6)
    plot_params$show_curves = FALSE
    plot_params$title = 'Scenario 1'
    model_params$plot_params = plot_params
    
    # Evaluation parameter list 
    eval_func_list = list()
    eval_func_list$mismatch = get_mismatches
    eval_func_list$vmeasure = get_v_measure

    model.result <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)
    true_m_not = model.result$`True basis coefficients`
    m_list = model.result$`M list`
    cluster_assignments = model.result$cluster_assignments
    MSE_list[[i]] <- get_MSE_1_2(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments)
  }
  list_sum <- Reduce("+", MSE_list)
  EMSE_k_tj <- list_sum / repeated_time
  return(EMSE_k_tj)
}

get_EMISE_1 <- function(repeated_time, T, n = 100){
  EMSE_k_tj <- get_EMSE_1(repeated_time)
  EMISE <- apply(EMSE_k_tj, 1, sum) 
  EMISE <- EMISE * T / n
  return(list(EMSE_k_t = EMSE_k_tj, EMISE = EMISE))
}

scenario.1.EMISE <- get_EMISE_1(50, T = pi/3, n = 100)

## For scenario 2--------------------------------------------------------------------------------------------------------------
# get MSE for each repetition
get_EMSE_2 <- function(repeated_time){
  MSE_list <- list()
  for (i in 1:repeated_time) {
    number_of_simulations = 1
    save_path = NULL
    
    # Data Parameters
    x = seq(from=0,to=pi/3, length = 100)
    K = 3
    curves_per_cluster = 50 
    true_cluster_assignments = rep(1:K,each = curves_per_cluster)
    seeds = i
    
    # Pack into data parameter list
    data_params = list()
    data_params$x = x
    data_params$K = K
    data_params$curves_per_cluster = curves_per_cluster
    data_params$true_cluster_assignments = true_cluster_assignments
    data_params$seeds = seeds 
    data_params$generate_data = Case_9
    
    # Model Parameters
    init = "km"
    nbasis = 6
    gamma_dist_config_matrix = matrix(0, 2, K)
    gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 25
    gamma_dist_config_matrix[2, ] = c(7.031, 7.031, 7.031) * 25
    d_not_vector = c(1/3, 1/3, 1/3)
    m_not_vector = matrix(c(0.45, 0.52, 0.69, 0.81, 0.66, 0.50,
                            0.73, 0.83, 1.07, 1.39, 1.51, 1.54,
                            0.62, 0.74, 1.03, 1.51, 1.84, 2.01), nrow = 3, byrow = T)
    v_not_vector = c(10, 10, 10, 10, 10, 10)
    convergence_threshold = 0.01
    max_iterations = 100
    verbose = FALSE
    draw = FALSE
    
    # Pack into model parameter list 
    model_params = list()
    model_params$model_func = get_funclustVI_cluster_assignments
    model_params$init = "km"
    model_params$nbasis = 6
    model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
    model_params$d_not_vector = d_not_vector
    model_params$m_not_vector = m_not_vector 
    model_params$v_not_vector = v_not_vector   
    model_params$convergence_threshold = convergence_threshold
    model_params$max_iterations = max_iterations 
    model_params$save_path = save_path
    model_params$verbose = verbose
    model_params$draw = draw
    plot_params = list()
    plot_params$xlim = NULL
    plot_params$ylim = c(1.5, 5)
    plot_params$show_curves = FALSE
    plot_params$title = 'Scenario 2'
    model_params$plot_params = plot_params
    
    # Evaluation parameter list 
    eval_func_list = list()
    eval_func_list$mismatch = get_mismatches
    eval_func_list$vmeasure = get_v_measure
    
    model.result <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)
    true_m_not = model.result$`True basis coefficients`
    m_list = model.result$`M list`
    cluster_assignments = model.result$cluster_assignments
    MSE_list[[i]] <- get_MSE_1_2(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments)
  }
  list_sum <- Reduce("+", MSE_list)
  EMSE_k_tj <- list_sum / repeated_time
  return(EMSE_k_tj)
}

get_EMISE_2 <- function(repeated_time, T, n = 100){
  EMSE_k_tj <- get_EMSE_2(repeated_time)
  EMISE <- apply(EMSE_k_tj, 1, sum) 
  EMISE <- EMISE * T / n
  return(list(EMSE_k_t = EMSE_k_tj, EMISE = EMISE))
}

scenario.2.EMISE <- get_EMISE_2(50, T = pi/3, n = 100)


## For scenario 3 and 4--------------------------------------------------------------------------------------------------------------
# get MSE for each repetition
get_MSE_3_4 <- function(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments){
  B <- get_B(x, nbasis)
  K <- nrow(true_m_not)
  MSE <- matrix(rep(0, K * length(x)), nrow = K)
  for (k in 1:K) {
    last_coeff <- c(m_list[[1]][1], m_list[[2]][1], m_list[[3]][1])
    j <- which.min((last_coeff-true_m_not[k, 1])^2)
    est_values <- as.vector(B %*% t(m_list[[j]]))
    true_values <- as.vector(B %*% true_m_not[k, ])
    MSE[k, ] <- (est_values - true_values)^2 
  }
  return(MSE)
}

# repetition
get_EMSE_3 <- function(repeated_time){
  MSE_list <- list()
  for (i in 1:repeated_time) {
    number_of_simulations = 1
    save_path = NULL
    
    # Data Parameters
    x <- seq(from=0, to=1, by=0.01)
    K = 3
    curves_per_cluster = 50 
    true_cluster_assignments = rep(1:K,each = curves_per_cluster)
    seeds = i
    
    # Pack into data parameter list
    data_params = list()
    data_params$x = x
    data_params$K = K
    data_params$curves_per_cluster = curves_per_cluster
    data_params$true_cluster_assignments = true_cluster_assignments
    data_params$seeds = seeds 
    data_params$generate_data = Case_1
    
    # Model Parameters
    init = "km"
    nbasis = 6
    gamma_dist_config_matrix = matrix(0, 2, K)
    gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 10
    gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 10
    d_not_vector = c(1/3, 1/3, 1/3)
    m_not_vector = matrix(c(1.5, 1, 1.8, 2.0, 1, 1.5,
                            2.8, 1.4, 1.8, 0.5, 1.5, 2.5,
                            0.4, 0.6, 2.4, 2.6, 0.1, 0.4), nrow = 3, byrow = T)
    v_not_vector = c(10, 10, 10, 10, 10, 10)
    convergence_threshold = 0.01
    max_iterations = 100
    verbose = FALSE
    draw = FALSE
    
    # Pack into model parameter list 
    model_params = list()
    model_params$model_func = get_funclustVI_cluster_assignments
    model_params$init = "km"
    model_params$nbasis = 6
    model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
    model_params$d_not_vector = d_not_vector
    model_params$m_not_vector = m_not_vector 
    model_params$v_not_vector = v_not_vector  
    model_params$convergence_threshold = convergence_threshold
    model_params$max_iterations = max_iterations 
    model_params$save_path = save_path
    model_params$verbose = verbose
    model_params$draw = draw
    plot_params = list()
    plot_params$xlim = NULL
    plot_params$ylim = c(-1, 4)
    plot_params$show_curves = FALSE
    plot_params$title = 'Scenario 3'
    model_params$plot_params = plot_params
    
    # Evaluation parameter list 
    eval_func_list = list()
    eval_func_list$mismatch = get_mismatches
    eval_func_list$vmeasure = get_v_measure
    
    model.result <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)
    true_m_not = model.result$`True basis coefficients`
    m_list = model.result$`M list`
    cluster_assignments = model.result$cluster_assignments
    MSE_list[[i]] <- get_MSE_3_4(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments)
  }
  list_sum <- Reduce("+", MSE_list)
  EMSE_k_tj <- list_sum / repeated_time
  return(EMSE_k_tj)
}

# for scenario 3
get_EMISE_3 <- function(repeated_time, T, n = 100){
  EMSE_k_tj <- get_EMSE_3(repeated_time)
  EMISE <- apply(EMSE_k_tj, 1, sum) 
  EMISE <- EMISE * T / n
  return(list(EMSE_k_t = EMSE_k_tj, EMISE = EMISE))
}

scenario.3.EMISE <- get_EMISE_3(50, T = 1, n = 100)


# for scenario 4
get_EMSE_4 <- function(repeated_time){
  MSE_list <- list()
  for (i in 1:repeated_time) {
    number_of_simulations = 1
    save_path = NULL
    
    # Data Parameters
    x <- seq(from=0, to=1, by=0.01)
    K = 3
    curves_per_cluster = 50 
    true_cluster_assignments = rep(1:K,each = curves_per_cluster)
    seeds = i
    
    # Pack into data parameter list
    data_params = list()
    data_params$x = x
    data_params$K = K
    data_params$curves_per_cluster = curves_per_cluster
    data_params$true_cluster_assignments = true_cluster_assignments
    data_params$seeds = seeds 
    data_params$generate_data = Case_2
    
    # Model Parameters
    init = "km"
    nbasis = 6
    gamma_dist_config_matrix = matrix(0, 2, K)
    gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 10
    gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 10
    d_not_vector = c(1/3, 1/3, 1/3)
    m_not_vector = matrix(c(1.5, 1, 1.6, 1.8, 1, 1.5,
                            1.8, 0.6, 0.4, 2.6, 2.8, 1.6,
                            1.2, 1.8, 2.2, 0.8, 0.6, 1.8), nrow = 3, byrow = T)
    v_not_vector = c(10, 10, 10, 10, 10, 10)
    convergence_threshold = 0.01
    max_iterations = 100
    verbose = FALSE
    draw = FALSE
    
    # Pack into model parameter list 
    model_params = list()
    model_params$model_func = get_funclustVI_cluster_assignments
    model_params$init = "km"
    model_params$nbasis = 6
    model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
    model_params$d_not_vector = d_not_vector
    model_params$m_not_vector = m_not_vector 
    model_params$v_not_vector = v_not_vector  
    model_params$convergence_threshold = convergence_threshold
    model_params$max_iterations = max_iterations 
    model_params$save_path = save_path
    model_params$verbose = verbose
    model_params$draw = draw
    plot_params = list()
    plot_params$xlim = NULL
    plot_params$ylim = c(-1, 4)
    plot_params$show_curves = FALSE
    plot_params$title = 'Scenario 4'
    model_params$plot_params = plot_params
    
    # Evaluation parameter list 
    eval_func_list = list()
    eval_func_list$mismatch = get_mismatches
    eval_func_list$vmeasure = get_v_measure
    
    model.result <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)
    true_m_not = model.result$`True basis coefficients`
    m_list = model.result$`M list`
    cluster_assignments = model.result$cluster_assignments
    MSE_list[[i]] <- get_MSE_3_4(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments)
  }
  list_sum <- Reduce("+", MSE_list)
  EMSE_k_tj <- list_sum / repeated_time
  return(EMSE_k_tj)
}

get_EMISE_4 <- function(repeated_time, T, n = 100){
  EMSE_k_tj <- get_EMSE_4(repeated_time)
  EMISE <- apply(EMSE_k_tj, 1, sum) 
  EMISE <- EMISE * T / n
  return(list(EMSE_k_t = EMSE_k_tj, EMISE = EMISE))
}

scenario.4.EMISE <- get_EMISE_4(50, T = 1, n = 100)

## For scenario 5--------------------------------------------------------------------------------------------------------------
# get MSE for each repetition
get_MSE_5 <- function(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments){
  B <- get_B(x, nbasis)
  K <- nrow(true_m_not)
  MSE <- matrix(rep(0, K * length(x)), nrow = K)
  for (k in 1:K) {
    last_coeff <- c(m_list[[1]][4], m_list[[2]][4], m_list[[3]][4])
    j <- which.min((last_coeff-true_m_not[k, 4])^2)
    est_values <- as.vector(B %*% t(m_list[[j]]))
    true_values <- as.vector(B %*% true_m_not[k, ])
    MSE[k, ] <- (est_values - true_values)^2 
  }
  return(MSE)
}

# repetition
get_EMSE_5 <- function(repeated_time){
  MSE_list <- list()
  for (i in 1:repeated_time) {
    # Initializationw
    number_of_simulations = 1
    save_path = NULL
    
    # Data Parameters
    x <- seq(from = 0, to = 24, length.out = 96)
    K = 3
    curves_per_cluster = 50 
    true_cluster_assignments = rep(1:K,each = curves_per_cluster)
    seeds = i
    
    # Pack into data parameter list
    data_params = list()
    data_params$x = x
    data_params$K = K
    data_params$curves_per_cluster = curves_per_cluster
    data_params$true_cluster_assignments = true_cluster_assignments
    data_params$seeds = seeds 
    data_params$generate_data = Case_5
    
    # Model Parameters
    init = "km"
    nbasis = 12
    gamma_dist_config_matrix = matrix(0, 2, K)
    gamma_dist_config_matrix[1, ] = c(347.22, 347.22, 347.22) 
    gamma_dist_config_matrix[2, ] = c(0.05, 0.05, 0.05) 
    d_not_vector = c(1/3, 1/3, 1/3)
    m_not_vector = matrix(c( 0.03, 0.07, -0.03, 0.19, 0.07, 0.05, 0.07, 0.03, 0.12, 0.05, 0.04, 0.04,
                             0.02, 0.01, 0.03, 0.17, -0.01, 0.03, 0.01, 0.03, 0.05, 0.01, 0.02, 0.02,
                             0.03, 0.03, 0.18, 0.02, 0.02, 0.02, 0.02, 0.06, 0.02, 0.02, 0.02, 0.02), nrow = 3, byrow = T)
    v_not_vector = rep(10, 12)
    convergence_threshold = 0.01
    max_iterations = 100
    verbose = FALSE
    draw = FALSE
    
    # Pack into model parameter list 
    model_params = list()
    model_params$model_func = get_funclustVI_cluster_assignments
    model_params$init = "km"
    model_params$nbasis = 12
    model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
    model_params$d_not_vector = d_not_vector
    model_params$m_not_vector = m_not_vector 
    model_params$v_not_vector = v_not_vector
    model_params$convergence_threshold = convergence_threshold
    model_params$max_iterations = max_iterations 
    model_params$save_path = save_path
    model_params$verbose = verbose
    model_params$draw = draw
    plot_params = list()
    plot_params$xlim = NULL
    plot_params$ylim =c(-0.05, 0.20)
    plot_params$show_curves = FALSE
    plot_params$title = 'Scenario 5'
    model_params$plot_params = plot_params
    
    # Evaluation parameter list 
    eval_func_list = list()
    eval_func_list$mismatch = get_mismatches
    eval_func_list$vmeasure = get_v_measure
    
    model.result <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)
    true_m_not = model.result$`True basis coefficients`
    m_list = model.result$`M list`
    cluster_assignments = model.result$cluster_assignments
    MSE_list[[i]] <- get_MSE_5(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments)
  }
  list_sum <- Reduce("+", MSE_list)
  EMSE_k_tj <- list_sum / repeated_time
  return(EMSE_k_tj)
}

get_EMISE_5 <- function(repeated_time, T, n = 96){
  EMSE_k_tj <- get_EMSE_5(repeated_time)
  EMISE <- apply(EMSE_k_tj, 1, sum) 
  EMISE <- EMISE * T / n
  return(list(EMSE_k_t = EMSE_k_tj, EMISE = EMISE))
}

scenario.5.EMISE <- get_EMISE_5(50, T = 24, n = 96)


## For scenario 6--------------------------------------------------------------------------------------------------------------
# get MSE for each repetition
get_MSE_6 <- function(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments){
  B <- get_B(x, nbasis)
  K <- nrow(true_m_not)
  MSE <- matrix(rep(0, K * length(x)), nrow = K)
  for (k in 1:K) {
    last_coeff <- c(m_list[[1]][4], m_list[[2]][4], m_list[[3]][4], m_list[[4]][4])
    j <- which.min((last_coeff-true_m_not[k, 4])^2)
    est_values <- as.vector(B %*% t(m_list[[j]]))
    true_values <- as.vector(B %*% true_m_not[k, ])
    MSE[k, ] <- (est_values - true_values)^2 
  }
  return(MSE)
}

# repetition
get_EMSE_6 <- function(repeated_time){
  MSE_list <- list()
  for (i in 1:repeated_time) {
    # Initializationw
    number_of_simulations = 1
    save_path = NULL
    
    # Data Parameters
    x = seq(from=0,to=pi/3, length = 100)
    K = 4
    curves_per_cluster = 50 
    true_cluster_assignments = rep(1:K,each = curves_per_cluster)
    seeds = i
    
    # Pack into data parameter list
    data_params = list()
    data_params$x = x
    data_params$K = K
    data_params$curves_per_cluster = curves_per_cluster
    data_params$true_cluster_assignments = true_cluster_assignments
    data_params$seeds = seeds 
    data_params$generate_data = Case_11
    
    # Model Parameters
    init = "km"
    nbasis = 6
    gamma_dist_config_matrix = matrix(0, 2, K)
    gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125, 78.125) * 25
    gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5, 12.5) * 25
    d_not_vector = c(1/4, 1/4, 1/4, 1/4)
    m_not_vector = matrix(c(0.06,   -0.34,   -1.13,   -0.54,    0.93,    1.66,
                            0.69,    0.18,   -0.78,    0.81,    2.44,    2.82,
                            0.64,    0.03,   -0.96,    1.43,    2.65,    2.62,
                            1.57,    0.83,   -0.08,    3.05,     3.4,    3.03), nrow = 4, byrow = T)
    v_not_vector = c(10, 10, 10, 10, 10, 10)
    convergence_threshold = 0.01
    max_iterations = 100
    verbose = FALSE
    draw = FALSE
    
    # Pack into model parameter list 
    model_params = list()
    model_params$model_func = get_funclustVI_cluster_assignments
    model_params$init = "km"
    model_params$nbasis = 6
    model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
    model_params$d_not_vector = d_not_vector
    model_params$m_not_vector = m_not_vector 
    model_params$v_not_vector = v_not_vector
    model_params$convergence_threshold = convergence_threshold
    model_params$max_iterations = max_iterations 
    model_params$save_path = save_path
    model_params$verbose = verbose
    model_params$draw = draw
    plot_params = list()
    plot_params$xlim = NULL
    plot_params$ylim = c(0, 6)
    plot_params$show_curves = FALSE
    plot_params$title = 'Scenario 6'
    model_params$plot_params = plot_params
    
    # Evaluation parameter list 
    eval_func_list = list()
    eval_func_list$mismatch = get_mismatches
    eval_func_list$vmeasure = get_v_measure
    
    model.result <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)
    true_m_not = model.result$`True basis coefficients`
    m_list = model.result$`M list`
    cluster_assignments = model.result$cluster_assignments
    MSE_list[[i]] <- get_MSE_6(x, nbasis, true_m_not, m_list, cluster_assignments, true_cluster_assignments)
  }
  list_sum <- Reduce("+", MSE_list)
  EMSE_k_tj <- list_sum / repeated_time
  return(EMSE_k_tj)
}

get_EMISE_6 <- function(repeated_time, T, n = 100){
  EMSE_k_tj <- get_EMSE_6(repeated_time)
  EMISE <- apply(EMSE_k_tj, 1, sum) 
  EMISE <- EMISE * T / n
  return(list(EMSE_k_t = EMSE_k_tj, EMISE = EMISE))
}

scenario.6.EMISE <- get_EMISE_6(50, T = pi/3, n = 100)

## Plots--------------------------------------------------------------------------------------------------------------
# plot Scenario 1
library(tidyverse)
scenario.1.EMSE.dt <- as.data.frame(t(scenario.1.EMISE$EMSE_k_t))
scenario.1.EMSE.dt$x <- seq(from=0,to=pi/3, length = 100)
colnames(scenario.1.EMSE.dt) <- c("Cluster 1", "Cluster 2", "Cluster 3", "x")
scenario.1.EMSE.dt %>% tidyr::gather("id", "EMSE", 1:3) %>% 
  ggplot(., aes(x, EMSE))+
  geom_line()+
  facet_wrap(~id) +
  ggtitle("The plot of EMSE versus the observed time point (x) in Scenario 1")
  
# plot Scenario 2
scenario.2.EMSE.dt <- as.data.frame(t(scenario.2.EMISE$EMSE_k_t))
scenario.2.EMSE.dt$x <- seq(from=0,to=pi/3, length = 100)
colnames(scenario.2.EMSE.dt) <- c("Cluster 1", "Cluster 2", "Cluster 3", "x")
scenario.2.EMSE.dt %>% tidyr::gather("id", "EMSE", 1:3) %>% 
  ggplot(., aes(x, EMSE))+
  geom_line()+
  facet_wrap(~id) +
  ggtitle("The plot of EMSE versus the observed time point (x) in Scenario 2")

# plot Scenario 3
scenario.3.EMSE.dt <- as.data.frame(t(scenario.3.EMISE$EMSE_k_t))
scenario.3.EMSE.dt$x <- seq(from=0,to=1, length = 101)
colnames(scenario.3.EMSE.dt) <- c("Cluster 1", "Cluster 2", "Cluster 3", "x")
scenario.3.EMSE.dt %>% tidyr::gather("id", "EMSE", 1:3) %>% 
  ggplot(., aes(x, EMSE))+
  geom_line()+
  facet_wrap(~id) +
  ggtitle("The plot of EMSE versus the observed time point (x) in Scenario 3")

# plot Scenario 4
scenario.4.EMSE.dt <- as.data.frame(t(scenario.4.EMISE$EMSE_k_t))
scenario.4.EMSE.dt$x <- seq(from=0,to=1, length = 101)
colnames(scenario.4.EMSE.dt) <- c("Cluster 1", "Cluster 2", "Cluster 3", "x")
scenario.4.EMSE.dt %>% tidyr::gather("id", "EMSE", 1:3) %>% 
  ggplot(., aes(x, EMSE))+
  geom_line()+
  facet_wrap(~id) +
  ggtitle("The plot of EMSE versus the observed time point (x) in Scenario 4")

# plot Scenario 5
scenario.5.EMSE.dt <- as.data.frame(t(scenario.5.EMISE$EMSE_k_t))
scenario.5.EMSE.dt$x <- seq(from=0,to=24, length = 96)
colnames(scenario.5.EMSE.dt) <- c("Cluster 1", "Cluster 2", "Cluster 3", "x")
scenario.5.EMSE.dt %>% tidyr::gather("id", "EMSE", 1:3) %>% 
  ggplot(., aes(x, EMSE))+
  geom_line()+
  facet_wrap(~id) +
  ggtitle("The plot of EMSE versus the observed time point (x) in Scenario 5")

# plot Scenario 6
scenario.6.EMSE.dt <- as.data.frame(t(scenario.6.EMISE$EMSE_k_t))
scenario.6.EMSE.dt$x <- seq(from=0,to=pi/3, length = 100)
colnames(scenario.6.EMSE.dt) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "x")
scenario.6.EMSE.dt %>% tidyr::gather("id", "EMSE", 1:4) %>% 
  ggplot(., aes(x, EMSE))+
  geom_line()+
  facet_wrap(~id) +
  ggtitle("The plot of EMSE versus the observed time point (x) in Scenario 6")

