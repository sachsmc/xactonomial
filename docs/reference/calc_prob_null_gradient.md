# Gradient of the multinomial likelihood sum

Gradient of the multinomial likelihood sum

## Usage

``` r
calc_prob_null_gradient(theta_cands, SSpacearr, II)
```

## Arguments

- theta_cands:

  A matrix with samples in the rows and the parameters in the columns

- SSpacearr:

  A matrix with the sample space for the given size of the problem

- II:

  logical vector of statistic for elements of sample space statistic
  being more extreme than the observed statistic

## Value

A matrix the same dimension as theta_cands

## Examples

``` r
calc_prob_null_gradient(t(c(.28, .32, .4)),
matrix(c(2, 2, 1, 1, 2, 2, 0, 3, 2), ncol = 3),
rep(TRUE, 3))
#>          [,1]      [,2]       [,3]
#> [1,] 0.199254 0.1102833 0.02679112

# numerically
testenv <- new.env()
testenv$SSpacearr <- matrix(c(2, 2, 1, 1, 2, 2, 0, 3, 2), ncol = 3)
testenv$thistheta <- c(.28, .32, .4)
numericDeriv(quote(sum(exp((.colSums(t(SSpacearr) * log(thistheta), m = 3, n = 3))))),
    theta = "thistheta", rho = testenv, central = TRUE)
#> [1] 0.03018932
#> attr(,"gradient")
#>          [,1]      [,2]       [,3]
#> [1,] 0.199254 0.1102833 0.02679112

```
