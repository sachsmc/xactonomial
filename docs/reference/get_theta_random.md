# Sample from the unit simplex in d dimensions

Sample from the unit simplex in d dimensions

## Usage

``` r
get_theta_random(d = 4, nsamp = 75)
```

## Arguments

- d:

  the dimension

- nsamp:

  the number of samples to take uniformly in the d space

## Value

The grid over Theta, the parameter space. A matrix with d columns and
nsamp rows

## Examples

``` r
get_theta_random(3, 10)
#>              [,1]       [,2]       [,3]
#>  [1,] 0.080750138 0.75358290 0.16566696
#>  [2,] 0.157208442 0.44355244 0.39923911
#>  [3,] 0.007399441 0.45899406 0.53360650
#>  [4,] 0.289767245 0.20801014 0.50222261
#>  [5,] 0.732881987 0.03963952 0.22747849
#>  [6,] 0.174940627 0.69966003 0.12539934
#>  [7,] 0.034241333 0.28614440 0.67961427
#>  [8,] 0.195669835 0.20665840 0.59767176
#>  [9,] 0.063661457 0.33987666 0.59646188
#> [10,] 0.388701313 0.58684652 0.02445216
```
