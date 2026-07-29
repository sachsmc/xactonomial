# Find a univariate root of the function f

This finds the value \\x \in \[a, b\]\\ such that \\f(x) = 0\\ using the
one-dimensional root finding ITP method (Interpolate Truncate Project).
Also see itp.

## Usage

``` r
itp_root(
  f,
  a,
  b,
  k1 = 0.1,
  k2 = 2,
  n0 = 1,
  eps = 0.005,
  maxiter = 100,
  fa = NULL,
  fb = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- f:

  The function to find the root of in terms of its first
  (one-dimensional) argument

- a:

  The lower limit

- b:

  The upper limit

- k1:

  A tuning parameter

- k2:

  Another tuning parameter

- n0:

  Another tuning parameter

- eps:

  Convergence tolerance

- maxiter:

  Maximum number of iterations

- fa:

  The value of f(a), if NULL then will be calculated

- fb:

  The value of f(b), if NULL then will be calculated

- verbose:

  Prints out information during iteration

- ...:

  Other arguments passed on to f

## Value

A numeric vector of length 1, the root at the last iteration

## References

I. F. D. Oliveira and R. H. C. Takahashi. 2020. An Enhancement of the
Bisection Method Average Performance Preserving Minmax Optimality. ACM
Trans. Math. Softw. 47, 1, Article 5 (March 2021), 24 pages.
https://doi.org/10.1145/3423597

## Examples

``` r
fpoly <- function(x) x^3 - x - 2 ## example from the ITP_method wikipedia entry
itp_root(fpoly, 1, 2, eps = .0001, verbose = TRUE)
#> iteration:  0 candidate:  1.433333 , closest value:  -0.4886296 
#> iteration:  1 candidate:  1.527131 , closest value:  0.03433833 
#> iteration:  2 candidate:  1.520093 , closest value:  -0.007641477 
#> iteration:  3 candidate:  1.521379 , closest value:  -4.253635e-06 
#> iteration:  4 candidate:  1.521383 , closest value:  1.964979e-05 
#> [1] 1.521381
#> attr(,"iter")
#> [1] 4
```
