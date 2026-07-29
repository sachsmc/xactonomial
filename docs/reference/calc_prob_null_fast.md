# Calculate probability for given parameters

Given a set of candidate parameter vectors, the enumerated sample space,
and a logical vector with the same number of elements of the sample
space, compute the probability for each element of the sample space and
take the sum.

## Usage

``` r
calc_prob_null_fast(theta_cands, SSpacearr, logC, II)
```

## Arguments

- theta_cands:

  A matrix with samples in the rows and the parameters in the columns

- SSpacearr:

  A matrix with the sample space for the given size of the problem

- logC:

  log multinomial coefficient for each element of the sample space

- II:

  logical vector of statistic for elements of sample space statistic
  being more extreme than the observed statistic

## Value

A numeric vector of probabilities

## Examples

``` r
sspace_3_5 <- matrix(sspace_multinom(3, 5), ncol = 3, byrow = TRUE)
theta_cands <- matrix(sample_unit_simplexn(3, 10), ncol = 3,byrow = TRUE)
calc_prob_null_fast(theta_cands, sspace_3_5,
apply(sspace_3_5, 1, log_multinom_coef, sumx = 5), II = 1:21 > 12)
#>  [1] 0.3851906695 0.2467206393 0.3338235041 0.0004533876 0.0076767952
#>  [6] 0.7663105384 0.7588514445 0.2469748988 0.4370544499 0.9424411110
# same as below but faster
calc_prob_null(theta_cands, sspace_3_5,
apply(sspace_3_5, 1, log_multinom_coef, sumx = 5), II = 1:21 > 12)
#>  [1] 0.3851906695 0.2467206393 0.3338235041 0.0004533876 0.0076767952
#>  [6] 0.7663105384 0.7588514445 0.2469748988 0.4370544499 0.9424411110
```
