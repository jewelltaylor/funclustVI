
<!-- README.md is generated from README.Rmd. Please edit that file -->

# funclustVI

The following package was created by John Jewell as a part of my
Undergraduate Dissertation at Western University under the supervision
of Professor Camila de Souza. The package serves to cluster functional
data using variational inference. More details in regards to the
functionality of the package is available in the examples section.

<!-- badges: start -->

<!-- badges: end -->

The goal of funclustVI is to …

## Installation

You can install the released version of funclustVI from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("funclustVI")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jewelltaylor/funclustVI")
```

## Example

This is an example which shows you how to use the package to generate
cluster assignments from functional data.

``` r
library(funclustVI)

# Model Arguments 
x = seq(from=0,to=pi/3, length = 100)
curves_per_cluster = 10
Y = Case_7(x, curves_per_cluster)
K = 3 
init = "km"
nbasis = 6
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 100
gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 100
convergence_threshold = 1
verbose = FALSE

model = funcslustVI(Y, K, nbasis, x, init, true_cluster_assignments, gamma_dist_config_matrix, convergence_threshold, verbose)

cluster_assignemnts = model$cluster_assignments

print(cluster_assignemnts)
#>  [1] 1 2 2 1 2 2 2 2 2 1 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1
```

This is an example which shows how to run simulations. Since draw =
True, a plot is generated showing the true vs estimated functions.

``` r
#Simulation Arguments 
case = "Case_7"
x = seq(from=0,to=pi/3, length = 100)
curves_per_cluster = 10
number_of_simulations = 1
Y = Case_7(x, curves_per_cluster)
K = 3 
initialization_method = "km"
nbasis = 6
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 100
gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 100
convergence_threshold = 1
save_path = NULL
verbose = FALSE
draw = TRUE 

simulation_results = run_simulations(case, x, K, curves_per_cluster, number_of_simulations, initialization_method, nbasis, true_cluster_assignments, gamma_dist_config_matrix, convergence_threshold, save_path, verbose, draw)
```

<img src="man/figures/README-simulation_plot-1.png" width="100%" />

    #> mm =  2 vm =  0.8410911 seed =  1 
    #> Average Mismatches:  2 
    #> Average V-measure:  0.8410911
