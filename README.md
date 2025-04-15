
# xactonomial

<!-- badges: start -->
<!-- badges: end -->

The goal of xactonomial is to use an exact (but computational and stochastic) method to compute a confidence interval and a function for calculation of p values in the k sample multinomial setting where interest is about a real-valued functional of the multinomial probabilities. 

## Installation

You can install the development version of xactonomial like so:

``` r
remotes::install_github("sachsmc/xactonomial")
```

but building the package from source requires a local rust environment. Instead, pre-built binaries can be found at https://sachsmc.r-universe.dev/xactonomial or used directly from R as
``` r
install.packages("xactonomial", repos = c('https://sachsmc.r-universe.dev', 'https://cloud.r-project.org'))
```


## Example

This is a basic example which shows you how to use the main function:

``` r
library(xactonomial)

psi_ba <- function(theta) {
   theta1 <- theta[1:4]
   theta2 <- theta[5:8]
   sum(sqrt(theta1 * theta2))
   }

data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
xactonomial(data, psi_ba, psi_limits = c(0, 1), maxit = 5, chunksize = 20)
```

