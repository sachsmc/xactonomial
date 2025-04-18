# Generated by extendr: Do not edit by hand

# nolint start

#
# This file was created with the following call:
#   .Call("wrap__make_xactonomial_wrappers", use_symbols = TRUE, package_name = "xactonomial")

#' @usage NULL
#' @useDynLib xactonomial, .registration = TRUE
NULL

#' Sample n times from the unit simplex in d dimensions
#' @param d the dimension
#' @param n the number of samples to take uniformly in the d space
#' @returns The grid over Theta, the parameter space. To be converted to a matrix with d columns and nsamp rows
#' @export
#' @examples
#' matrix(sample_unit_simplexn(3, 10), ncol = 3, byrow = TRUE)
sample_unit_simplexn <- function(d, n) .Call(wrap__sample_unit_simplexn, d, n)

#' Calculate multinomial probabilities
#' @param sar The unrolled matrix containing the portion of the sample space to sum over
#' @param logt The vector of candidate theta values, as sampled from the null space
#' @param logc The vector of log multinomial coefficients see \link{log_multinom_coef}
#' @param d The total dimension, sum(d_j)
#' @param n The sample size
#' @param nt The number of candidate theta values
#' @returns A vector of probabilities
#' @export
#' @examples
#' sspace_3_5 <- sspace_multinom(3, 5)
#' calc_multinom_probs(sspace_3_5, sample_unit_simplexn(3, 10),
#'   apply(matrix(sspace_3_5, ncol = 3, byrow = TRUE), 1, log_multinom_coef, sumx = 5), 3, 5, 10)
#'
calc_multinom_probs <- function(sar, logt, logc, d, n, nt) .Call(wrap__calc_multinom_probs, sar, logt, logc, d, n, nt)

#' Enumerate the multinomial sample space
#' @param d The dimension
#' @param n The sample size
#' @returns A vector enumerating the sample space, to be converted to a matrix
#' with d columns and choose(n + d - 1, d - 1) rows
#' @export
#' @examples
#' matrix(sspace_multinom(3, 5), ncol = 3, byrow = TRUE)
sspace_multinom <- function(d, n) .Call(wrap__sspace_multinom, d, n)


# nolint end
