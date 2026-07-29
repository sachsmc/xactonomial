# Enumerate the multinomial sample space

Enumerate the multinomial sample space

## Usage

``` r
sspace_multinom(d, n)
```

## Arguments

- d:

  The dimension

- n:

  The sample size

## Value

A vector enumerating the sample space, to be converted to a matrix with
d columns and choose(n + d - 1, d - 1) rows

## Examples

``` r
matrix(sspace_multinom(3, 5), ncol = 3, byrow = TRUE)
#>       [,1] [,2] [,3]
#>  [1,]    5    0    0
#>  [2,]    4    1    0
#>  [3,]    3    2    0
#>  [4,]    2    3    0
#>  [5,]    1    4    0
#>  [6,]    0    5    0
#>  [7,]    4    0    1
#>  [8,]    3    1    1
#>  [9,]    2    2    1
#> [10,]    1    3    1
#> [11,]    0    4    1
#> [12,]    3    0    2
#> [13,]    2    1    2
#> [14,]    1    2    2
#> [15,]    0    3    2
#> [16,]    2    0    3
#> [17,]    1    1    3
#> [18,]    0    2    3
#> [19,]    1    0    4
#> [20,]    0    1    4
#> [21,]    0    0    5
```
