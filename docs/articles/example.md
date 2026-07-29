# Basic Examples

``` r

library(xactonomial)
```

We consider the following problem. Let $`T_j`$ be distributed
$`\mbox{Multinomial}(n_j, d_j, \theta_j)`$ for $`j = 1, \ldots, k`$,
where $`n_j`$ is the number of trials, $`d_j`$ the number of mutually
exclusive categories, and $`\theta_j`$ the vector of $`d_j`$
probabilities. Denote $`\boldsymbol{T} = (T_1, \ldots, T_k)`$ and
$`\boldsymbol{\theta} = (\theta_1, \ldots, \theta_k)`$. Suppose one is
interested in the parameter
$`\psi \equiv \tau(\boldsymbol{\theta}) \in \Psi`$ where $`\Psi`$ is a
subset of the real line $`\mathbb{R}`$. Suppose we observe a sample
$`\boldsymbol{X} = (X_1, \ldots, X_k)`$ which is a vector of counts
corresponding to the random variable $`\boldsymbol{T}`$. Let
$`G(\boldsymbol{X})`$ denote a real-valued statistic that defines an
ordering of the sample space. A sensible default statistic in this
setting is to apply the function defining the parameter of interest to
the sample proportions.

If the function $`\tau()`$ is complex, or the variance of
$`\psi(\hat{\boldsymbol{\theta}})`$ is otherwise difficult to calculate,
then one may consider using the nonparametric bootstrap for inference
(Efron 1979). However, the bootstrap may perform poorly in small
samples, e.g., not having the correct coverage or type I error (Bickel
and Freedman 1981). Instead we propose a computational approach for
computing p values and confidence intervals in this setting.

Here is a simple illustration of the problem. Let
``` math
\psi = \tau(\theta_1, \theta_2) = \sum \sqrt(\theta_1, \theta_2), 
```
which is similar to the Bhattacharyya coefficient (Bhattacharyya 1946).
Suppose we have two samples $`T_1`$ and $`T_2`$ each from independent
multinomials with dimension 4. The following R code allows us to sample
the data and compute $`\psi`$.

``` r

true_theta <- c(.45, .15, .3, .1, .05, .15, .4, .4)

sample_data <- function(n) {
  
  T1 <- rmultinom(1, n[1], prob = true_theta[1:4])
  T2 <- rmultinom(1, n[2], prob = true_theta[5:8])
  
  list(T1 = c(T1), T2 = c(T2))
  
}

tau_bc <- function(theta) {
  
  theta1 <- theta[1:4]
  theta2 <- theta[5:8]
  sum(sqrt(theta1 * theta2))
  
}

tau_bc(true_theta)
#> [1] 0.8464102
```

We can do the bootstrap, but it actually performs quite poorly in this
situation. We are going to define the same psi function but vectorized
this time. This improves the speed of xactonomial.

``` r


tau_bc_v <- function(theta) {
  
  theta1 <- theta[,1:4, drop = FALSE]
  theta2 <- theta[,5:8, drop = FALSE]
  rowSums(sqrt(theta1 * theta2))
  
}
```

``` r

set.seed(510)
data <- sample_data(n = c(10, 12))
data
#> $T1
#> [1] 9 0 1 0
#> 
#> $T2
#> [1] 1 1 4 6
results <- xactonomial(data, tau_bc_v, psi_limits = c(0,1), conf_int = TRUE, psi0 = .5,
                          maxit = 50, chunksize = 100, ga_restart_every = 10)
results
#> 
#>  Monte Carlo multinomial test
#> 
#> data:  data
#> p-value = 0.7103
#> alternative hypothesis: true psi0 is not equal to 0.5
#> 95 percent confidence interval:
#>  0.3192390 0.9284745
#> sample estimates:
#>  tau_bc_v 
#> 0.4564355
```

## Maximum of parameters

In this example, $`\tau(\theta) = \max(\theta)`$. When two of the
parameters are equal, the function is nondifferentiable, so bootstrap
and other asymptotic methods will fail while our method still works.

``` r

tau_max <- function(pp) {
  
 max(pp)
  
}

true_tau <- tau_max(c(.4, .4, .2))

sample_data2 <- function(n) {
  
  list(rmultinom(1, n, prob = c(.4, .4, .2))[, 1])
  
}

set.seed(421)
data <- sample_data2(n = c(60))
data
#> [[1]]
#> [1] 27 25  8
results <- xactonomial(data, tau_max, psi_limits = c(1/3,1), psi0 = .55,
                          maxit = 100, chunksize = 1000, itp_eps = .01)
#> Warning in xactonomial(data, tau_max, psi_limits = c(1/3, 1), psi0 = 0.55, : No
#> root found for lower confidence limit, using lower boundary.

results
#> 
#>  Monte Carlo multinomial test
#> 
#> data:  data
#> p-value = 0.154
#> alternative hypothesis: true psi0 is not equal to 0.55
#> 95 percent confidence interval:
#>  0.3333333 0.5883333
#> sample estimates:
#> tau_max 
#>    0.45
```

We can get rid of the warning by providing the argument `p_value_limits`
which should be the vector of 2 p-values, the first being for the test
of the null at the lower boundary of psi and the second at the upper
boundary. If those p-values are greater than $`\alpha / 2`$, then no
root exists and the confidence limit should instead be the boundary.

First we calculate the exact p-value for the test of `psi0 <= 1/3` by
providing `theta_null_points`. The upper p-value is clearly less than
$`\alpha / 2`$.

``` r

p.ll <- xactonomial(data, tau_max, psi_limits = c(1/3,1), psi0 = 1/3, conf_int = FALSE,
   itp_eps = .01, theta_null_points = t(rep(1/3, 3)), 
   alternative = "greater")$p.value

xactonomial(data, tau_max, psi_limits = c(1/3,1), psi0 = 1, conf_int = FALSE,
   itp_eps = .01, theta_null_points = rbind(c(1, 0, 0), c(0,1,0), c(0,0,1)), 
   alternative = "less")$p.value
#> [1] 0
```

Then we provide the vector of those 2 p-values to the function. They do
not need to be precise, the important thing is when one of them is
greater than $`\alpha / 2`$.

``` r

xactonomial(data, tau_max, psi_limits = c(1/3,1), p_value_limits = c(p.ll, 1e-8),
                          maxit = 100, chunksize = 1000, itp_eps = .01)
#> 
#>  Monte Carlo multinomial test
#> 
#> data:  data
#> p-value = NA
#> alternative hypothesis: two.sided
#> 95 percent confidence interval:
#>  0.3333333 0.5883333
#> sample estimates:
#> tau_max 
#>    0.45
```

Now you can see we get exactly 1/3 for the lower confidence limit.

## References

Bhattacharyya, A. 1946. “On a Measure of Divergence Between Two
Multinomial Populations.” *Sankhyā: The Indian Journal of Statistics
(1933-1960)* 7 (4): 401–6. <http://www.jstor.org/stable/25047882>.

Bickel, Peter J, and David A Freedman. 1981. “Some Asymptotic Theory for
the Bootstrap.” *The Annals of Statistics* 9 (6): 1196–217.

Efron, B. 1979. “Bootstrap Methods: Another Look at the Jackknife.” *The
Annals of Statistics* 7 (1): 1–26.
