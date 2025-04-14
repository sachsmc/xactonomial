#' Enumerate the sample space of a multinomial
#'
#' We have \eqn{d} mutually exclusive outcomes and \eqn{n} independent trials.
#' This function enumerates all possible vectors of length \eqn{d} of counts of
#' each outcome for \eqn{n} trials, i.e., the sample space. The result is output
#' as a matrix with \eqn{d} columns where each row represents a possible
#' observation.
#'
#' @param d Dimension
#' @param n Size
#' @returns A matrix with d columns
#' @export
#' @examples
#' d4s <- sspace_multinom(4, 8)
#' stopifnot(abs(sum(apply(d4s, 1, dmultinom, prob = rep(.25, 4))) - 1) < 1e-12)
#'

sspace_multinom1 <- function(d, n) {

  if(d == 2) {

    cbind(0:n, n:0)

  } else {

    res <- NULL
    for(i in 0:n) {

      res <- rbind(res, cbind(i, sspace_multinom(d - 1, n - i)))

    }
    res

  }

}


sspace_multinom2 <- function(d, n) {

  ncomb <- choose(d + n - 1, d - 1)
  binss <- matrix(rep(c(n, rep(0, d - 1)), ncomb), nrow = ncomb, byrow = TRUE)
  bins <- c(n, rep(0, d - 1))
  i <- 2
  repeat{
    if(bins[d] == n) break
    if(bins[1] > 0) {
      bins[1] <- bins[1] - 1
      bins[2] <- bins[2] + 1
      #i <- i + 1
    } else {
      nz <- 2
      while(bins[nz] == 0) {
        nz <- nz + 1
      }
        bins[1] <- bins[nz] - 1
        bins[nz + 1] <- bins[nz + 1] + 1
        bins[nz] <- 0
        #i <- i + 1
      }
    binss[i, ] <- bins
    i <- i + 1
  }

  binss

}



#' Get a matrix of indices for all possible combinations of vectors of lengths
#'
#' This is basically the same as \link[base]{expand.grid}, but faster for integers
#'
#' @param lengths A vector with the lengths of each index to expand
#' @returns A matrix with length(lengths) columns and prod(lengths) rows
#' @export
#' @examples
#' expand_index(c(2, 3, 4))
#' ## the same as
#' expand.grid(1:2, 1:3, 1:4)
#'
expand_index <- function(lengths) {

  do.call(expand.grid, lapply(lengths, seq_len)) |>
    as.matrix()

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




#' Sample from the unit simplex in d dimensions
#' @param d the dimension
#' @param nsamp the number of samples to take uniformly in the d space
#' @return The grid over Theta, the parameter space. A matrix with d columns and nsamp rows
#' @export
#' @examples
#' get_theta_random(3, 10)
#'

get_theta_random <- function(d = 4, nsamp = 75) {

  # x1 <- matrix(runif(nsamp * (d - 1)), ncol = d - 1)
  # g2 <- unique(t(apply(x1, 1, \(x) {
  #   diff(sort(c(0, x, 1)))
  # })))

  matrix(sample_unit_simplexn(d, nsamp), nrow = nsamp, ncol = d, byrow = TRUE)

}


#' Sample uniformly from d_k simplexes
#' @param d_k vector of vector lengths
#' @param nsamp number of samples to take
#' @returns A matrix with sum(d_k) columns and nsamp rows
#' @export
#' @examples
#' runif_dk_vects(c(3, 4, 2), 10)
runif_dk_vects <- function(d_k, nsamp, ...){

  do.call("cbind", lapply(d_k, \(i) get_theta_random(i, nsamp)))

}

#' Sample from a Dirichlet distribution for each of d_k vectors
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

    #sample_dirichlet(nsamp, ceiling(alpha[[i]]))
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

  (a + b) / 2


}




#' Calculate combinations of multinomials
#'
#' @param X Matrix 1
#' @param Y Matrix 2
#' @returns A list of arrays
#'
#' @export
combinate <- function(X, Y) {
  Yi <- rep(1:nrow(Y), rep.int(nrow(X), nrow(Y)))
  Xi <- rep(1:nrow(X), times = ceiling(length(Yi)/nrow(X)))
  newX <- X[Xi,]
  newY <- Y[Yi,]

  sumX <- sum(newX[1,])
  sumY <- sum(newY[1,])

  logCX <- lfactorial(sumX) - rowSums(apply(newX, 2, lfactorial))
  logCY <- lfactorial(sumY) - rowSums(apply(newY, 2, lfactorial))

  list(Sspace = cbind(newX, newY), Sprobs = cbind(newX / sumX, newY / sumY), logC = logCX + logCY)

}


#' Like combinate but add to existing
#'
#' @param X A list as returned by combinate
#' @param Y Matrix 2
#' @returns A list of arrays
#'
#' @export
combinate2 <- function(X, Y) {
  Yi <- rep(1:nrow(Y), rep.int(nrow(X$Sspace), nrow(Y)))
  Xi <- rep(1:nrow(X$Sspace), times = ceiling(length(Yi)/nrow(X$Sspace)))
  newX <- X$Sspace[Xi,]
  newY <- Y[Yi,]

  sumY <- sum(newY[1,])
  logCY <- lfactorial(sumY) - rowSums(apply(newY, 2, lfactorial))

  list(Sspace = cbind(newX, newY), Sprobs = cbind(X$Sprobs[Xi,], newY / sumY), logC = X$logC[Xi] + logCY)

}
