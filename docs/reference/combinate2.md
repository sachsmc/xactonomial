# Like [combinate](https://sachsmc.github.io/xactonomial/reference/combinate.md) but adds on to previous call

Like
[combinate](https://sachsmc.github.io/xactonomial/reference/combinate.md)
but adds on to previous call

## Usage

``` r
combinate2(X, Y)
```

## Arguments

- X:

  A list containing the elements Sspace (matrix), and logC (vector), the
  result of a call to
  [combinate](https://sachsmc.github.io/xactonomial/reference/combinate.md)

- Y:

  Matrix 2

## Value

A list containing Sspace, the sample space (vectors of counts), and
logC, a vector of the log multinomial coefficients.

## Examples

``` r
slist_2_3 <- combinate(matrix(sspace_multinom(2, 5), ncol = 2, byrow = TRUE),
   matrix(sspace_multinom(3, 6), ncol = 3, byrow = TRUE))

sl_2_3_4 <- combinate2(slist_2_3, matrix(sspace_multinom(4, 3), ncol = 4, byrow = TRUE))
```
