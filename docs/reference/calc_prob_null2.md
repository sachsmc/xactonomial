# Calculate probability for given parameters

Given a set of candidate parameter vectors, check if the null \\\psi
\leq \psi_0\\ is satisfied, and if so, compute the probability for each
element of the sample space

## Usage

``` r
calc_prob_null2(
  theta_cands,
  psi,
  psi0,
  minus1,
  SSpacearr,
  logC,
  II,
  psi_v = FALSE
)
```

## Arguments

- theta_cands:

  A matrix with samples in the rows and the parameters in the columns

- psi:

  The function of interest mapping parameters to the real line

- psi0:

  The null boundary for testing psi \<= psi0

- minus1:

  Either plus or minus 1

- SSpacearr:

  A matrix with the sample space for the given size of the problem

- logC:

  log multinomial coefficient for each element of the sample space

- II:

  logical vector of sample space psi being more extreme than the
  observed psi

- psi_v:

  Is psi vectorized by row?

## Value

A numeric vector of probabilities
