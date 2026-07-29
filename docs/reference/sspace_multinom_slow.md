# Enumerate the sample space of a multinomial

We have \\d\\ mutually exclusive outcomes and \\n\\ independent trials.
This function enumerates all possible vectors of length \\d\\ of counts
of each outcome for \\n\\ trials, i.e., the sample space. The result is
output as a matrix with \\d\\ columns where each row represents a
possible observation. See
[sspace_multinom](https://sachsmc.github.io/xactonomial/reference/sspace_multinom.md)
for a faster implementation using Rust.

## Usage

``` r
sspace_multinom_slow(d, n)
```

## Arguments

- d:

  Dimension

- n:

  Size

## Value

A matrix with d columns

## Examples

``` r
d4s <- sspace_multinom_slow(4, 8)
stopifnot(abs(sum(apply(d4s, 1, dmultinom, prob = rep(.25, 4))) - 1) < 1e-12)
```
