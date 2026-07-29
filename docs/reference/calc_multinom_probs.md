# Calculate multinomial probabilities

Calculate multinomial probabilities

## Usage

``` r
calc_multinom_probs(sar, logt, logc, d, n, nt)
```

## Arguments

- sar:

  The unrolled matrix containing the portion of the sample space to sum
  over

- logt:

  The vector of candidate theta values, as sampled from the null space

- logc:

  The vector of log multinomial coefficients see
  [log_multinom_coef](https://sachsmc.github.io/xactonomial/reference/log_multinom_coef.md)

- d:

  The total dimension, sum(d_j)

- n:

  The sample size

- nt:

  The number of candidate theta values

## Value

A vector of probabilities

## Examples

``` r
sspace_3_5 <- sspace_multinom(3, 5)
calc_multinom_probs(sspace_3_5, sample_unit_simplexn(3, 10),
  apply(matrix(sspace_3_5, ncol = 3, byrow = TRUE), 1, log_multinom_coef, sumx = 5), 3, 5, 10)
#>  [1] 296.84007 142.02929 106.60159 108.70758 295.14131 304.71594  69.98053
#>  [8]  84.68910  86.55943 357.01053
```
