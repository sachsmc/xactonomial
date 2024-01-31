#' Enumerate the sample space of a multinomial of dimension d with size n
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
  lgamma(sumx+1) - sum( lgamma(x+1) )
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


#' Null probability for parameters
#'
#' Given a set of candidate parameter vectors, check if the null is satisfied, and if so, compute the probability for each element of the sample space
#'
#' @param theta_cands A matrix with samples in the rows and the parameters in the columns
#' @param psi The function of interest mapping parameters to the real line
#' @param psi0 The null boundary for testing psi <= psi0
#' @param minus1 Either plus or minus 1
#' @param SSpace A list with the relevant components
#' @param II logical vector of sample space psi being more extreme than the observed
#' @param logC log multinomial coefficient for each element of the sample space
#'
#' @returns A numeric vector
#'
#' @export
#'

calc_prob_null <- function(theta_cands, psi, psi0, minus1, SSspacearr, II, logC) {

  k <- dim(SSpacearr)[2]
  res <- combn(1:nrow(theta_cands), k, FUN = \(i) {

    thisthetav <- lapply(i, \(ii) theta_cands[ii,])
    thispsi <- do.call(psi, thisthetav)
    thetavect <- do.call(rbind, thisthetav)

    if(minus1 * thispsi <= minus1 * psi0) {

      sum(exp((apply(SSpacearr, 3, \(SSpacemat) sum(diag(log(thetavect) %*% SSpacemat))) + logC)[II]))
      #sum(exp((outer(c(S1 %*% log(theta1)), c(S2 %*% log(theta2)), "+") + logC)[II]))

    } else {
      NA
    }
  })
  res[!is.na(res)]

}


#' Get a matrix of indices for all possible combinations of vectors of lengths
#'
#' @param lengths A vector with the lengths of each index to expand
#' @returns A matrix with length(lengths) columns and prod(lengths) rows
#' @export
expand_index <- function(lengths) {

  orep <- prod(lengths)
  cdex <- matrix(NA, nrow = orep, ncol = length(lengths))

  for(i in 1:length(lengths)) {
    cdex[,i] <- rep.int(rep.int(seq_len(lengths[i]), rep.int(1,
                                         lengths[i])), orep / lengths[i])
  }

  cdex

}
