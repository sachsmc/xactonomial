# Compute a p value for the test of psi \<= psi0 and/or psi \>= psi0

Compute a p value for the test of psi \<= psi0 and/or psi \>= psi0

## Usage

``` r
pvalue_psi0(
  psi0,
  f_param,
  psi_hat,
  psi_obs,
  alternative = "two.sided",
  maxit,
  chunksize,
  p_target,
  SSpacearr,
  logC,
  d_k,
  f_is_vectorized = FALSE,
  theta_sampler = runif_dk_vects,
  ga = FALSE,
  ga_gfactor = 1,
  ga_lrate = 0.01,
  ga_restart_every = 10,
  warn = TRUE
)
```

## Arguments

- psi0:

  The null hypothesis value for the parameter being tested.

- f_param:

  Function that takes in parameters and outputs a real valued number for
  each parameter. Can be vectorized rowwise for a matrix or not.

- psi_hat:

  The vector of psi values at each element of the sample space

- psi_obs:

  The observed estimate at the given data

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of "two.sided" (default), "greater" or "less"

- maxit:

  Maximum number of iterations of the Monte Carlo procedure

- chunksize:

  The number of samples to take from the parameter space at each
  iteration

- p_target:

  If a p-value is found that is greater than p_target, terminate the
  algorithm early.

- SSpacearr:

  The sample space matrix

- logC:

  The log multinomial coefficients for each row of the sample space

- d_k:

  The vector of dimensions

- f_is_vectorized:

  Is f_param vectorized by row?

- theta_sampler:

  Function to take samples from the \\Theta\\ parameter space. Default
  is
  [runif_dk_vects](https://sachsmc.github.io/xactonomial/reference/runif_dk_vects.md).

- ga:

  Logical, if TRUE, uses gradient ascent.

- ga_gfactor:

  Concentration parameter scale in the gradient ascent algorithm. A
  number or "adapt"

- ga_lrate:

  The gradient ascent learning rate

- ga_restart_every:

  Restart the gradient ascent after this number of iterations at a
  sample from

- warn:

  If TRUE, will give a warning if no samples from the null space are
  found

## Value

A vector with two p-values, one for the lower, and one for the greater

## Examples

``` r

sspace_3_5 <- matrix(sspace_multinom(3, 5), ncol = 3, byrow = TRUE)
f_max <- function(theta) max(theta)
logC <- apply(sspace_3_5, 1, log_multinom_coef, sumx = 5)
psi_hat <- apply(sspace_3_5, 1, \(x) f_max(x / sum(x)))
pvalue_psi0(.3, f_max, psi_hat, .4, maxit = 10, chunksize = 100,
 p_target = 1, SSpacearr = sspace_3_5, logC = logC, d_k = 3, warn = FALSE)
#>      null       alt 
#>        NA 0.3702166 
#> attr(,"p.sequence")
#> attr(,"p.sequence")$p.null
#>  [1] 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
#> 
#> attr(,"p.sequence")$p.alt
#>  [1] 0.3650079 0.3681825 0.3681825 0.3681825 0.3681825 0.3702166 0.3702166
#>  [8] 0.3702166 0.3702166 0.3702166
#> 
```
