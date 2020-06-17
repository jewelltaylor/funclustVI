
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
install.packages("fda")
library(fda)
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jewelltaylor/funclustVI")
```

## Example

This is a basic example which shows you how to use the package to
generate cluster assignments from functional data/

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
#>  [1] 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 3 1 1 1
```

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
