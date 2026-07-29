# Improved inference for a real-valued function of multinomial parameters

We consider the k sample multinomial problem where we observe k vectors
(possibly of different lengths), each representing an independent sample
from a multinomial. For a given function \\\tau(\theta)\\ which takes in
the concatenated vector of multinomial probabilities \\\theta\\ and
outputs a real number, we are interested in computing a p-value for a
test of \\\tau(\theta) = \psi \geq \psi_0\\, and constructing a
confidence interval for \\\psi\\.

## Usage

``` r
xactonomial(
  data,
  f_param,
  statistic = NULL,
  psi0 = NULL,
  alternative = c("two.sided", "less", "greater"),
  psi_limits,
  theta_null_points = NULL,
  p_target = 1,
  conf_int = TRUE,
  conf_level = 0.95,
  itp_maxit = 10,
  itp_eps = 0.005,
  p_value_limits = NULL,
  maxit = 50,
  chunksize = 500,
  theta_sampler = runif_dk_vects,
  ga = TRUE,
  ga_gfactor = "adapt",
  ga_lrate = 0.01,
  ga_restart_every = 10,
  seed = 503
)
```

## Arguments

- data:

  A list with k elements representing the vectors of counts of a
  k-sample multinomial

- f_param:

  Function that takes in parameters and outputs psi, a real valued
  number for each parameter. Can be vectorized rowwise for a matrix or
  not.

- statistic:

  Function that takes in a matrix with data vectors in the rows, and
  outputs a vector with the number of rows in the matrix. If NULL, will
  be inferred from `f_param` by plugging in the empirical proportions.

- psi0:

  The null hypothesis value for the parameter being tested.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of "two.sided" (default), "greater" or "less"

- psi_limits:

  A vector of length 2 giving the lower and upper limits of the range of
  \\\tau(\theta)\\

- theta_null_points:

  An optional matrix where each row is a theta value that gives
  `f_param(theta) = psi0`. If this is supplied and psi0 = one of
  `psi_limits`, then a truly exact p-value will be calculated.

- p_target:

  If a p-value is found that is greater than p_target, terminate the
  algorithm early.

- conf_int:

  If TRUE, calculates a confidence interval by inverting the p-value
  function

- conf_level:

  A number between 0 and 1, the confidence level.

- itp_maxit:

  Maximum iterations to use in the ITP algorithm. Only relevant if
  conf_int = TRUE.

- itp_eps:

  Epsilon value to use for the ITP algorithm. Only relevant if conf_int
  = TRUE.

- p_value_limits:

  A vector of length 2 giving lower bounds on the one-sided p-values
  corresponding to psi0 = psi_limits, with alternative = "less" and
  "greater", respectively. Only relevant if conf_int = TRUE. See
  examples.

- maxit:

  Maximum number of iterations of the Monte Carlo procedure

- chunksize:

  The number of samples to take from the parameter space at each
  iteration

- theta_sampler:

  Function to take samples from the \\Theta\\ parameter space. Default
  is
  [runif_dk_vects](https://sachsmc.github.io/xactonomial/reference/runif_dk_vects.md).
  Must be a function of two parameters `d_k` a vector of dimensions, and
  `chunksize` the number of samples to take, and return a matrix with
  `sum(d_k)` columns and `chunksize` rows. See examples.

- ga:

  Logical, if TRUE, uses gradient ascent.

- ga_gfactor:

  Concentration parameter scale in the gradient ascent algorithm. A
  number or "adapt"

- ga_lrate:

  The gradient ascent learning rate

- ga_restart_every:

  Restart the gradient ascent after this number of iterations at a
  sample from `theta_sampler`

- seed:

  Seed for the random number generator. Can be set to NULL in which case
  no seed is set.

## Value

An object of class "htest", which is a list with the following elements:

- estimate:

  The value of the statistic at the observed data

- p.value:

  The p value

- conf.int:

  The upper and lower confidence limits

- null.value:

  The null hypothesis value provided by the user

- alternative:

  The type of test

- method:

  A description of the method

- data.name:

  The name of the data object provided by the user

- p.sequence:

  A list with two elements, p.null and p.alt containing the vector of p
  values at each iteration for the less than null and the greater than
  null. Used for assessing convergence.

## Details

Let \\T_j\\ be distributed
\\\mbox{Multinomial}\_{d_j}(\boldsymbol{\theta}\_j, n_j)\\ for \\j = 1,
\ldots, k\\ and denote \\\boldsymbol{T} = (T_1, \ldots, T_k)\\ and
\\\boldsymbol{\theta} = (\theta_1, \ldots, \theta_k)\\. The subscript
\\d_j\\ denotes the dimension of the multinomial. Suppose one is
interested in the parameter \\\psi = \tau(\boldsymbol{\theta}) \in \Psi
\subseteq \mathbb{R}\\. Given a sample of size \\n\\ from
\\\boldsymbol{T}\\, say \\\boldsymbol{X} = (X_1, \ldots, X_k)\\, which
is a vector of counts obtained by concatenating the k independent count
vectors, let \\G(\boldsymbol{X})\\ denote a real-valued statistic that
defines the ordering of the sample space. The default choice of the
statistic is to estimate \\\boldsymbol{\theta}\\ with the sample
proportions and plug them into \\\tau(\boldsymbol{\theta})\\. This
function calculates a p value for a test of the null hypothesis \\H_0:
\psi \neq \psi_0\\ for the two sided case, \\H_0: \psi \leq \psi_0\\ for
the case `alternative = "greater"`, and \\H_0: \psi \geq \psi_0\\ for
the case `alternative = "less"`. We make no assumptions and do not rely
on large sample approximations. It also optionally constructs a \\1 -
\alpha\\ percent confidence interval for \\\psi\\. The computation is
somewhat involved so it is best for small sample sizes. The calculation
is done by sampling a large number of points from the null parameter
space \\\Theta_0\\, then computing multinomial probabilities under those
values for the range of the sample space where the statistic is as or
more extreme than the observed statistic given data. It is basically the
definition of a p-value implemented with Monte Carlo methods. Some
options for speeding up the calculation are available.

## Specifying the function \\\tau(\cdot)\\

This function is the `f_param` argument and should be a function that
either: 1) takes a vector of length sum(d_j) (the total number of bins)
and outputs a single number, or 2) takes a matrix with number of columns
equal to sum(d_j), and arbitrary number of rows and outputs a vector
with length equal to the number of rows. In other words, psi can be not
vectorized or it can be vectorized by row. Writing it so that it is
vectorized can speed up the calculation. See examples.

## Boundary issues

It is required to provide `psi_limits`, a vector of length 2 giving the
smallest and largest possible values that the function psi can take,
e.g., `c(0, 1)`. If the null hypothesis value psi0 is at one of the
limits, it is often the case that sampling from the null parameter space
is impossible because it is a set of measure 0. While it may have
measure 0, it is not empty, and will contain a finite set of points.
Thus you should provide the argument `theta_null_points` which is a
matrix where the rows contain the finite set (sometimes 1) of points
\\\theta\\ such that \\\tau(\theta) = \psi_0\\. There is also an
argument called `p_value_limits` that can be used to improve performance
of confidence intervals around the boundary. This should be a vector of
length 2 with the p-value for a test of psi_0 \<= psi_limits\[1\] and
the p-value for a test of psi_0 \>= psi_limits\[2\]. See examples.

## Optimization options

For p-value calculation, you can provide a parameter p_target, so that
the sampling algorithm terminates when a p-value is found that exceeds
p_target. The algorithm begins by sampling uniformly from the unit
simplices defining the parameter space, but alternatives can be
specified in `theta_sampler`. By default gradient ascent (`ga = TRUE`)
is performed during the p-value maximization procedure, and `ga_gfactor`
and `ga_lrate` control options for the gradient ascent. At each
iteration, the gradient of the multinomial probability at the current
maximum theta is computed, and a step is taken to
`theta + lrate * gradient`. Then for the next iteration, a set of
`chunksize` samples are drawn from a Dirichlet distribution with
parameter `ga_gfactor * (theta + ga_lrate * gradient)`. If
`ga_gfactor = "adapt"` then it is set to `1 / max(theta)` at each
iteration. The ITP algorithm
[itp_root](https://sachsmc.github.io/xactonomial/reference/itp_root.md)
is used to find roots of the p-value function as a function of the psi0
value to get confidence intervals. The maximum number of iterations and
epsilon can be controlled via `itp_maxit, itp_eps`.

## References

Sachs, M.C., Gabriel, E.E. and Fay, M.P., 2024. Exact confidence
intervals for functions of parameters in the k-sample multinomial
problem. arXiv preprint arXiv:2406.19141.

## Examples

``` r
tau_ba <- function(theta) {
  theta1 <- theta[1:4]
  theta2 <- theta[5:8]
  sum(sqrt(theta1 * theta2))
  }
data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
xactonomial(data, tau_ba, psi_limits = c(0, 1), psi0 = .5,
  conf_int = FALSE, maxit = 15, chunksize = 200)
#> 
#>  Monte Carlo multinomial test
#> 
#> data:  data
#> p-value = 0.02401
#> alternative hypothesis: true psi0 is not equal to 0.5
#> 95 percent confidence interval:
#>  NA NA
#> sample estimates:
#>    tau_ba 
#> 0.7995291 
#> 

# vectorized by row
tau_ba_v <- function(theta) {
theta1 <- theta[,1:4, drop = FALSE]
theta2 <- theta[,5:8, drop = FALSE]
rowSums(sqrt(theta1 * theta2))
}
data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
xactonomial(data, tau_ba_v, psi_limits = c(0, 1), psi0 = .5,
 conf_int = FALSE, maxit = 10, chunksize = 200)
#> 
#>  Monte Carlo multinomial test
#> 
#> data:  data
#> p-value = 0.02401
#> alternative hypothesis: true psi0 is not equal to 0.5
#> 95 percent confidence interval:
#>  NA NA
#> sample estimates:
#>  tau_ba_v 
#> 0.7995291 
#> 

 # example of using theta_null_points
 # psi = 1/3 occurs when all probs = 1/3
 tau_max <- function(pp) {
   max(pp)
 }

data <- list(c(13, 24, 13))

xactonomial(data, tau_max, psi_limits = c(1 / 3, 1), psi0 = 1/ 3,
  conf_int = FALSE, theta_null_points = t(c(1/3, 1/3, 1/3)))
#> 
#>  Exact multinomial test given a point null
#> 
#> data:  data
#> p-value = 0.1331
#> alternative hypothesis: true psi0 is not equal to 0.3333333
#> 95 percent confidence interval:
#>  NA NA
#> sample estimates:
#> tau_max 
#>    0.48 
#> 

## in this case using p_value_limits improves confidence interval performance

 xactonomial(data, tau_max, psi_limits = c(1 / 3, 1), psi0 = 1/ 3,
  conf_int = TRUE, theta_null_points = t(c(1/3, 1/3, 1/3)),
  p_value_limits = c(.1, 1e-8))
#> 
#>  Exact multinomial test given a point null
#> 
#> data:  data
#> p-value = 0.1331
#> alternative hypothesis: true psi0 is not equal to 0.3333333
#> 95 percent confidence interval:
#>  0.3333333 0.6258333
#> sample estimates:
#> tau_max 
#>    0.48 
#> 

## specifying theta_sampler

dirich_sampler <- function(d_k, chunksize){
 rdirich_dk_vects(chunksize, list(1:4 + 1))
}

xactonomial(list(1:4),tau_max,
            psi_limits = c(0.25,1), psi0 = .5, conf_int = FALSE,
            theta_sampler = dirich_sampler)
#> 
#>  Monte Carlo multinomial test
#> 
#> data:  list(1:4)
#> p-value = 0.6586
#> alternative hypothesis: true psi0 is not equal to 0.5
#> 95 percent confidence interval:
#>  NA NA
#> sample estimates:
#> tau_max 
#>     0.4 
#> 
```
