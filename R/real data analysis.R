### real data analysis: Growth-------------------------------------------------------------------------------------------
library(fda)
boy.ht <- t(as.matrix(growth$hgtm))
girl.ht <- t(as.matrix(growth$hgtf))
ht.data <- rbind(boy.ht, girl.ht)
colnames(ht.data) <- NULL
rownames(ht.data) <- NULL

# Data Parameters
x = seq(1:31)
Y = ht.data 
K = 2
curves_per_cluster = c(39, 54)
true_cluster_assignments = rep(1:K, curves_per_cluster)

# Model Parameters 
init = "km"
nbasis = 10 
convergence_threshold = 0.01
max_iterations = 100
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(100, 100) * 20
gamma_dist_config_matrix[2, ] = c(5, 5) * 20
d_not_vector = c(1/3, 2/3)
m_not_vector = matrix(c(70,   82,   85,  122,  141,   148,  177,  180,  181,   181,
                        63,   78,   83,  118,  135,   140,  150,  158,  158,   158), 
                      nrow = 2, byrow = T)
v_not_vector = rep(10, 10)
verbose = FALSE
draw = TRUE
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(66, 196)
plot_params$show_curves = FALSE
plot_params$title = "Growth"

# Fit the model 
growth.clust = funcslustVI(x, Y, K, true_cluster_assignments, 
                    init, nbasis, convergence_threshold, max_iterations, 
                    gamma_dist_config_matrix, d_not_vector, m_not_vector, v_not_vector,
                    verbose, draw, plot_params)

Mismatch(growth.clust$cluster_assignments, true_cluster_assignments, 2)
sabre::vmeasure(growth.clust$cluster_assignments, true_cluster_assignments)$v_measure

# repetitions
growth.model <- function(seed){
  set.seed(seed)
  model = funcslustVI(x, Y, K, true_cluster_assignments, 
                      init, nbasis, convergence_threshold, max_iterations, 
                      gamma_dist_config_matrix, d_not_vector, m_not_vector, v_not_vector, 
                      verbose, draw, plot_params)
  
  mis.mat <- Mismatch(model$cluster_assignments, true_cluster_assignments, 2)
  v.measure <- sabre::vmeasure(model$cluster_assignments, true_cluster_assignments)$v_measure
  return(c(mis.mat, v.measure))
}


growth.analysis <- function(repetition.times){
  start.time <- Sys.time()
  growth.metrics <- sapply(1:repetition.times, growth.model)
  mis.mean <- mean(growth.metrics[1, ])
  v.mean <- mean(growth.metrics[2, ])
  mis.sd <- sd(growth.metrics[1, ])
  v.sd <- sd(growth.metrics[2, ])
  end.time <- Sys.time()
  time.length <- end.time - start.time
  return(list(mismatch.mean = mis.mean,
              vmeasure.mean = v.mean,
              mismatch.sd = mis.sd,
              vmeasure.sd = v.sd,
              runing.time = time.length))
}
growth.performance <- growth.analysis(50)

# plot of Growth
cluster_number = true_cluster_assignments[1]
col = 2 + cluster_number
plot(x, Y[1, ], col=col, type="l", ylim=c(66, 196), xlim=c(1, 31), main="Growth", ylab="f(x)", xlab="x")
for(i in 2:NROW(Y)) {
  cluster_number = true_cluster_assignments[i]
  col = 2 + cluster_number
  lines(x, Y[i, ], type="l", col=col)
}

# k-means repetitions for comparison
k.means.growth.model <- function(seed){
  set.seed(seed)
  model = kmeans(Y, 2)
  mis.mat <- Mismatch(model$cluster, true_cluster_assignments, 2)
  v.measure <- sabre::vmeasure(model$cluster, true_cluster_assignments)$v_measure
  return(c(mis.mat, v.measure))
}

k.means.growth.analysis <- function(repetition.times){
  start.time <- Sys.time()
  growth.metrics <- sapply(1:repetition.times, k.means.growth.model)
  mis.mean <- mean(growth.metrics[1, ])
  v.mean <- mean(growth.metrics[2, ])
  mis.sd <- sd(growth.metrics[1, ])
  v.sd <- sd(growth.metrics[2, ])
  end.time <- Sys.time()
  time.length <- end.time - start.time
  return(list(mismatch.mean = mis.mean,
              vmeasure.mean = v.mean,
              mismatch.sd = mis.sd,
              vmeasure.sd = v.sd,
              runing.time = time.length))
}
k.means.growth.performance <- k.means.growth.analysis(50)



### real data analysis: Canadian weather-------------------------------------------------------------------------------------------
wea.data <- t(as.matrix(CanadianWeather$dailyAv[,,1]))
nrow(wea.data)
ncol(wea.data)
class(CanadianWeather$region)
table(CanadianWeather$region)
colnames(wea.data) <- NULL
rownames(wea.data) <- NULL

# Data Parameters
x = seq(1:365)
Y = wea.data 
K = 4
curves_per_cluster = c(15, 12, 5, 3)
true_cluster_assignments = rep(1:K, curves_per_cluster)

# Model Parameters 
init = "km"
nbasis = 6
convergence_threshold = 0.001
max_iterations = 1000
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(100, 100, 100, 100) * 10
gamma_dist_config_matrix[2, ] = c(80, 80, 80, 80) * 10
d_not_vector = c(1/2, 1/3, 1/10, 1/15)
m_not_vector = matrix(c(-3,   -9,    3,   20,   -8,   -3,
                        -14,  -16,   10,   20,   -6,  -14,
                        -10,   -7,    9,   22,   -7,  -10,
                        -24,  -28,   -3,   15,  -20,  -24), 
                      nrow = 4, byrow = T)
v_not_vector = rep(10, 6)
verbose = FALSE
draw = TRUE
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(-35, 23)
plot_params$show_curves = FALSE
plot_params$title = "Canadian weather"

# Fit the model 
set.seed(2)
model = funcslustVI(x, Y, K, true_cluster_assignments, 
                    init, nbasis, convergence_threshold, max_iterations, 
                    gamma_dist_config_matrix, d_not_vector, m_not_vector, v_not_vector,
                    verbose, draw, plot_params)

Mismatch(model$cluster_assignments, true_cluster_assignments, 4)
sabre::vmeasure(model$cluster_assignments, true_cluster_assignments)$v_measure

# repetitions
wea.model <- function(seed){
  wea.data <- t(as.matrix(CanadianWeather$dailyAv[,,1]))
  colnames(wea.data) <- NULL
  rownames(wea.data) <- NULL
  
  #Data Parameters
  x = seq(1:365)
  Y = wea.data 
  K = 4
  curves_per_cluster = c(15, 12, 5, 3)
  true_cluster_assignments = rep(1:K, curves_per_cluster)
  #Model Parameters 
  init = "km"
  nbasis = 6
  convergence_threshold = 0.001
  max_iterations = 1000
  gamma_dist_config_matrix = matrix(0, 2, K)
  gamma_dist_config_matrix[1, ] = c(100, 100, 100, 100) * 10
  gamma_dist_config_matrix[2, ] = c(80, 80, 80, 80) * 10
  d_not_vector = c(1/2, 1/3, 1/10, 1/15)
  m_not_vector = matrix(c(-3,   -9,    3,   20,   -8,   -3,
                          -14,  -16,   10,   20,   -6,  -14,
                          -10,   -7,    9,   22,   -7,  -10,
                          -24,  -28,   -3,   15,  -20,  -24), 
                        nrow = 4, byrow = T)
  v_not_vector = rep(10, 6)
  verbose = FALSE
  draw = FALSE
  plot_params = list()
  plot_params$xlim = NULL
  plot_params$ylim = c(-35, 23)
  plot_params$show_curves = FALSE
  plot_params$title = NULL
  
  #Fit the model 
  set.seed(seed)
  model = funcslustVI(x, Y, K, true_cluster_assignments, 
                      init, nbasis, convergence_threshold, max_iterations, 
                      gamma_dist_config_matrix, d_not_vector, m_not_vector, v_not_vector,
                      verbose, draw, plot_params)
  
  mis.mat = Mismatch(model$cluster_assignments, true_cluster_assignments, 4)
  v.measure = sabre::vmeasure(model$cluster_assignments, true_cluster_assignments)$v_measure
  return(c(mis.mat, v.measure))
}

wea.analysis <- function(repetition.times){
  start.time <- Sys.time()
  wea.metrics <- sapply(1:(repetition.times), wea.model)
  mis.mean <- mean(wea.metrics[1, ])
  v.mean <- mean(wea.metrics[2, ])
  mis.sd <- sd(wea.metrics[1, ])
  v.sd <- sd(wea.metrics[2, ])
  end.time <- Sys.time()
  time.length <- end.time - start.time
  return(list(mismatch.mean = mis.mean,
              vmeasure.mean = v.mean,
              mismatch.sd = mis.sd,
              vmeasure.sd = v.sd,
              runing.time = time.length))
}
wea.performance <- wea.analysis(10)

# k-means repetitions for comparison
k.means.wea.model <- function(seed){
  set.seed(seed)
  model = kmeans(Y, 5)
  mis.mat <- Mismatch(model$cluster, true_cluster_assignments, 5)
  v.measure <- sabre::vmeasure(model$cluster, true_cluster_assignments)$v_measure
  return(c(mis.mat, v.measure))
}

k.means.wea.analysis <- function(repetition.times){
  start.time <- Sys.time()
  wea.metrics <- sapply(1:(repetition.times), k.means.wea.model)
  mis.mean <- mean(wea.metrics[1, ])
  v.mean <- mean(wea.metrics[2, ])
  mis.sd <- sd(wea.metrics[1, ])
  v.sd <- sd(wea.metrics[2, ])
  end.time <- Sys.time()
  time.length <- end.time - start.time
  return(list(mismatch.mean = mis.mean,
              vmeasure.mean = v.mean,
              mismatch.sd = mis.sd,
              vmeasure.sd = v.sd,
              runing.time = time.length))
}
k.means.wea.performance <- k.means.wea.analysis(10)
k.means.wea.performance

# plot of Weather
cluster_number = true_cluster_assignments[1]
col = 2 + cluster_number
plot(x, Y[1, ], col=col, type="l", ylim=c(-10, 20), xlim=c(1, 365), main="Canadian weather", ylab="f(x)", xlab="x")
for(i in 2:NROW(Y)) {
  cluster_number = true_cluster_assignments[i]
  col = 2 + cluster_number
  lines(x, Y[i, ], type="l", col=col)
}
