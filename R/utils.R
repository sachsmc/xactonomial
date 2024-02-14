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

sspace_multinom <- function(d, n) {

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



#' Get a matrix of indices for all possible combinations of vectors of lengths
#'
#' @param lengths A vector with the lengths of each index to expand
#' @returns A matrix with length(lengths) columns and prod(lengths) rows
#' @export
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

get_theta_random <- function(d = 4, nsamp = 75) {

  x1 <- matrix(runif(nsamp * (d - 1)), ncol = d - 1)
  g2 <- unique(t(apply(x1, 1, \(x) {
    diff(sort(c(0, x, 1)))
  })))

  g2


}



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





smart_sample_theta <- function(psi, psi0, psi_limits, d_k, chunksize = 200, lower = TRUE) {

  psinull <- if(lower) {
    runif(chunksize, psi_limits[1], psi0)
  } else {
    runif(chunksize, psi0, psi_limits[2])
  }

  f_theta <- \(theta, psi_j) {

    theta <- plogis(theta)
    start <- 1
    for(d in d_k) {
      theta[start:(start + d - 1)] <- theta[start:(start + d - 1)] / sum(theta[start:(start + d - 1)])
      start <- start + d
    }
    (psi(theta) - psi_j)^2

  }

  #randstart <- do.call(cbind, lapply(d_k, \(i) qlogis(get_theta_random(i, chunksize))))
  randstart <- matrix(rnorm(chunksize * sum(d_k), sd = 2), nrow = chunksize, ncol = sum(d_k))
  out_thetas <- matrix(NA, nrow = chunksize, ncol = ncol(randstart))
  for(i in 1:chunksize) {
    res <- optim(randstart[i,], f_theta, psi_j = psinull[i])
    theta_i <- res$par
    theta_i <- plogis(theta_i)
    start <- 1
    for(d in d_k) {
      theta_i[start:(start + d - 1)] <- theta_i[start:(start + d - 1)] / sum(theta_i[start:(start + d - 1)])
      start <- start + d
    }

    out_thetas[i, ] <- theta_i
  }

  #out_thetas
  chkpsi <- apply(out_thetas, 1, psi)
  if(lower) out_thetas[chkpsi <= psi0,] else out_thetas[chkpsi >= psi0, ]

}



