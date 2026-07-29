# Get a sumsample of the multinomial sample space

Get a sumsample of the multinomial sample space

## Usage

``` r
sspace_multinom_sample(d, n, k)
```

## Arguments

- d:

  The dimension

- n:

  The sample size

- k:

  The number of elements to keep

## Value

A vector with a subsample from the sample space, to be converted to a
matrix with d columns and k \<= choose(n + d - 1, d - 1) rows

## Examples

``` r
matrix(sspace_multinom_sample(3, 5, 10), ncol = 3, byrow = TRUE)
#>       [,1] [,2] [,3]
#>  [1,]    5    0    0
#>  [2,]    4    1    0
#>  [3,]    1    4    0
#>  [4,]    0    5    0
#>  [5,]    4    0    1
#>  [6,]    0    4    1
#>  [7,]    3    0    2
#>  [8,]    2    1    2
#>  [9,]    0    3    2
#> [10,]    0    2    3
#> [11,]    0    0    5
```
