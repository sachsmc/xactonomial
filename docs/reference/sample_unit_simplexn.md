# Sample n times from the unit simplex in d dimensions

Sample n times from the unit simplex in d dimensions

## Usage

``` r
sample_unit_simplexn(d, n)
```

## Arguments

- d:

  the dimension

- n:

  the number of samples to take uniformly in the d space

## Value

The grid over Theta, the parameter space. To be converted to a matrix
with d columns and nsamp rows

## Examples

``` r
matrix(sample_unit_simplexn(3, 10), ncol = 3, byrow = TRUE)
#>            [,1]      [,2]        [,3]
#>  [1,] 0.0312627 0.4731823 0.495555002
#>  [2,] 0.2865548 0.5979900 0.115455134
#>  [3,] 0.3199876 0.4730253 0.206987082
#>  [4,] 0.1287485 0.6915184 0.179733058
#>  [5,] 0.2400753 0.3965726 0.363352079
#>  [6,] 0.4336132 0.3507270 0.215659802
#>  [7,] 0.4848650 0.1705140 0.344620999
#>  [8,] 0.2615073 0.7300833 0.008409427
#>  [9,] 0.2430331 0.1227204 0.634246439
#> [10,] 0.1317241 0.6727937 0.195482257
```
