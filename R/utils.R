#' Enumerate the sample space of a multinomial
#'
#' We have \eqn{d} mutually exclusive outcomes and \eqn{n} independent trials.
#' This function enumerates all possible vectors of length \eqn{d} of counts of
#' each outcome for \eqn{n} trials, i.e., the sample space. The result is output
#' as a matrix with \eqn{d} columns where each row represents a possible
#' observation. See \link{sspace_multinom} for a faster implementation using Rust.
#'
#' @param d Dimension
#' @param n Size
#' @returns A matrix with d columns
#' @export
#' @examples
#' d4s <- sspace_multinom_slow(4, 8)
#' stopifnot(abs(sum(apply(d4s, 1, dmultinom, prob = rep(.25, 4))) - 1) < 1e-12)
#'

sspace_multinom_slow <- function(d, n) {

  if(d == 2) {

    cbind(0:n, n:0)

  } else {

    res <- NULL
    for(i in 0:n) {

      res <- rbind(res, cbind(i, sspace_multinom_slow(d - 1, n - i)))

    }
    res

  }

}


#' Calculate log of multinomial coefficient
#' @param x Vector of observed counts in each cell
#' @param size Total count
#' @returns The log multinomial coefficient
#' @export
#' @examples
#' #' @examples
#' S0 <- sspace_multinom(4, 6)
#' S1 <- sspace_multinom(4, 7)
#' logC0<- apply(S0,1,log_multinom_coef,sumx=6)
#' logC1<- apply(S1,1,log_multinom_coef,sumx=7)
#' logC<- outer(logC0,logC1,'+')
#'
log_multinom_coef<- function(x,sumx){
  lfactorial(sumx) - sum( lfactorial(x) )
}


#' Sample uniformly and independently from d_k simplices
#' @param d_k vector of vector lengths
#' @param nsamp number of samples to take
#' @returns A matrix with sum(d_k) columns and nsamp rows
#' @export
#' @examples
#' runif_dk_vects(c(3, 4, 2), 10)
runif_dk_vects <- function(d_k, nsamp, ...){

  do.call("cbind", lapply(d_k, \(i) matrix(sample_unit_simplexn(i, nsamp), nrow = nsamp, byrow = TRUE)))

}

#' Sample independently from Dirichlet distributions for each of d_k vectors
#' @param nsamp number of samples to take
#' @param alpha List of vectors of concentration parameters
#' @returns A matrix with sum(d_k) columns and nsamp rows
#' @export
#' @examples
#' rdirich_dk_vects(10, list(rep(1, 3), rep(1, 4), rep(1, 2)))
#'
rdirich_dk_vects <- function(nsamp, alpha) {

  d_k <- sapply(alpha, length)
  do.call("cbind", lapply(1:length(d_k), \(i) {

    gsamps <- rgamma(nsamp * length(alpha[[i]]), shape = rep(alpha[[i]], each = nsamp))
    yis <- matrix(gsamps, nrow = nsamp)
    yis / rowSums(yis)

  }))

}


#' Find the root of the function f
#'
#' This finds the value \eqn{x \in [a, b]} such that \eqn{f(x) = 0} using the one-dimensional root finding ITP method (Interpolate Truncate Project). Also see \link[itp]{itp}.
#'
#' @param f The function to find the root of in terms of its first (one-dimensional) argument
#' @param a The lower limit
#' @param b The upper limit
#' @param k1 A tuning parameter
#' @param k2 Another tuning parameter
#' @param n0 Another tuning parameter
#' @param eps Convergence tolerance
#' @param maxiter Maximum number of iterations
#' @param fa The value of f(a), if NULL then will be calculated
#' @param fb The value of f(b), if NULL then will be calculated
#' @param verbose Prints out information during iteration
#' @param ... Other arguments passed on to f
#'
#' @returns A numeric vector of length 1, the root at the last iteration
#'
#' @references I. F. D. Oliveira and R. H. C. Takahashi. 2020. An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality. ACM Trans. Math. Softw. 47, 1, Article 5 (March 2021), 24 pages. https://doi.org/10.1145/3423597
#'
#' @export
#' @examples
#' fpoly <- function(x) x^3 - x - 2 ## example from the ITP_method wikipedia entry
#' itp_root(fpoly, 1, 2, eps = .0001, verbose = TRUE)
#'
itp_root <- function(f, a, b, k1 = .1, k2 = 2, n0 = 1,
                     eps = .005, maxiter = 100, fa = NULL, fb = NULL,
                     verbose = FALSE, ...) {

  if(is.null(fa) | is.null(fb)) {
    fa <- f(a, ...)
    fb <- f(b, ...)
  }
  inc <- sign(fb)

  n12 <- log((b - a) / (2 * eps), base = 2)
  nmax <- n12 + n0
  j <- 0

  for_rk <- 2 ^ (n0 - 1 + log2(eps) + ceiling(log2(b - a) - log2(eps)))

  repeat {
    x12 <- (a + b) / 2
    xf <- (fb * a - fa * b) / (fb - fa)


    delta <- k1 * (b - a)^k2
    sigma <- sign(x12 - xf)

    if(delta <= abs(x12 - xf)) {
      xt <- xf + sigma * delta
    } else {
      xt <- x12
    }

    r <- for_rk - (b - a) / 2


    if(abs(xt - x12) <= r) {
      xitp <- xt
    } else {
      xitp <- x12 - sigma * r
    }

    yitp <- f(xitp, ...)
    if(yitp * inc > 0) {
      b <- xitp
      fb <- yitp
    } else if(yitp * inc < 0) {
      a <- xitp
      fa <- yitp
    } else {
      a <- b <- xitp
    }

    if(verbose) {
      cat("iteration: ", j, "candidate: ", xitp, ", closest value: ", yitp, "\n")
    }

    if((b - a) < 2 * eps | j >= maxiter) break

    for_rk <- for_rk * .5
    j <- j + 1

  }

  res <- (a + b) / 2
  attr(res, "iter") <- j
  res


}




#' Calculate combinations of multinomial vectors
#'
#' Given X and Y, both matrices where the rows are counts of multinomial trials,
#' produce all combinations rowwise, and calculate the log multinomial
#' coefficients for the combination.
#' @param X Matrix 1
#' @param Y Matrix 2
#' @returns A list containing Sspace, the sample space (vectors of counts), and
#'  logC, a vector of the log multinomial coefficients.
#'
#' @export
#' @examples
#' slist_2_3 <- combinate(matrix(sspace_multinom(2, 5), ncol = 2, byrow = TRUE),
#'    matrix(sspace_multinom(3, 6), ncol = 3, byrow = TRUE))
#'
combinate <- function(X, Y) {
  Yi <- rep(1:nrow(Y), rep.int(nrow(X), nrow(Y)))
  Xi <- rep(1:nrow(X), times = ceiling(length(Yi)/nrow(X)))
  newX <- X[Xi,]
  newY <- Y[Yi,]

  sumX <- sum(newX[1,])
  sumY <- sum(newY[1,])

  logCX <- lfactorial(sumX) - rowSums(apply(newX, 2, lfactorial))
  logCY <- lfactorial(sumY) - rowSums(apply(newY, 2, lfactorial))

  list(Sspace = cbind(newX, newY), logC = logCX + logCY)

}


#' Like \link{combinate} but adds on to previous call
#'
#' @param X A list containing the elements Sspace (matrix), and logC (vector)
#' @param Y Matrix 2
#' @returns A list containing Sspace, the sample space (vectors of counts), and
#'  logC, a vector of the log multinomial coefficients.
#'
#' @export
#' @examples
#' slist_2_3 <- combinate(matrix(sspace_multinom(2, 5), ncol = 2, byrow = TRUE),
#'    matrix(sspace_multinom(3, 6), ncol = 3, byrow = TRUE))
#'
#' sl_2_3_4 <- combinate2(slist_2_3, matrix(sspace_multinom(4, 3), ncol = 4, byrow = TRUE))
combinate2 <- function(X, Y) {
  Yi <- rep(1:nrow(Y), rep.int(nrow(X$Sspace), nrow(Y)))
  Xi <- rep(1:nrow(X$Sspace), times = ceiling(length(Yi)/nrow(X$Sspace)))
  newX <- X$Sspace[Xi,]
  newY <- Y[Yi,]

  sumY <- sum(newY[1,])
  logCY <- lfactorial(sumY) - rowSums(apply(newY, 2, lfactorial))

  list(Sspace = cbind(newX, newY), logC = X$logC[Xi] + logCY)

}
