#' Enumerate the sample space of a multinomial of dimension d with size n
#'
#' @param d Dimension
#' @param n Size
#' @returns A matrix with d columns
#'
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

get_theta_random <- function(d = 4, nsamp = 75) {

  x1 <- matrix(runif(nsamp * (d - 1)), ncol = d - 1)
  g2 <- unique(t(apply(x1, 1, \(x) {
    diff(sort(c(0, x, 1)))
  })))

  g2


}
