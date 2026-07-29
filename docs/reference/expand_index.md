# Get a matrix of indices for all possible combinations of vectors of lengths

This is basically the same as
[expand.grid](https://rdrr.io/r/base/expand.grid.html), but faster for
integers

## Usage

``` r
expand_index(lengths)
```

## Arguments

- lengths:

  A vector with the lengths of each index to expand

## Value

A matrix with length(lengths) columns and prod(lengths) rows

## Examples

``` r
expand_index(c(2, 3, 4))
#>       Var1 Var2 Var3
#>  [1,]    1    1    1
#>  [2,]    2    1    1
#>  [3,]    1    2    1
#>  [4,]    2    2    1
#>  [5,]    1    3    1
#>  [6,]    2    3    1
#>  [7,]    1    1    2
#>  [8,]    2    1    2
#>  [9,]    1    2    2
#> [10,]    2    2    2
#> [11,]    1    3    2
#> [12,]    2    3    2
#> [13,]    1    1    3
#> [14,]    2    1    3
#> [15,]    1    2    3
#> [16,]    2    2    3
#> [17,]    1    3    3
#> [18,]    2    3    3
#> [19,]    1    1    4
#> [20,]    2    1    4
#> [21,]    1    2    4
#> [22,]    2    2    4
#> [23,]    1    3    4
#> [24,]    2    3    4
## the same as
expand.grid(1:2, 1:3, 1:4)
#>    Var1 Var2 Var3
#> 1     1    1    1
#> 2     2    1    1
#> 3     1    2    1
#> 4     2    2    1
#> 5     1    3    1
#> 6     2    3    1
#> 7     1    1    2
#> 8     2    1    2
#> 9     1    2    2
#> 10    2    2    2
#> 11    1    3    2
#> 12    2    3    2
#> 13    1    1    3
#> 14    2    1    3
#> 15    1    2    3
#> 16    2    2    3
#> 17    1    3    3
#> 18    2    3    3
#> 19    1    1    4
#> 20    2    1    4
#> 21    1    2    4
#> 22    2    2    4
#> 23    1    3    4
#> 24    2    3    4
```
