# Calculate log of multinomial coefficient

Calculate log of multinomial coefficient

## Usage

``` r
log_multinom_coef(x, sumx)
```

## Arguments

- x:

  Vector of observed counts in each cell

- sumx:

  Total count

## Value

The vector of log multinomial coefficients

## Examples

``` r
S0 <- matrix(sspace_multinom(4, 6), ncol = 4, byrow = TRUE)
logC0<- apply(S0,1,log_multinom_coef,sumx=6)
```
