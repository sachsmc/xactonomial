#' Calculate probability for given parameters
#'
#' Given a set of candidate parameter vectors, check if the null \eqn{\psi \leq
#' \psi_0} is satisfied, and if so, compute the probability for each element of
#' the sample space
#'
#' @param theta_cands A matrix with samples in the rows and the parameters in
#'   the columns
#' @param psi The function of interest mapping parameters to the real line
#' @param psi0 The null boundary for testing psi <= psi0
#' @param minus1 Either plus or minus 1
#' @param SSpacearr A matrix with the sample space for the given size of the
#'   problem
#' @param II logical vector of sample space psi being more extreme than the
#'   observed psi
#' @param logC log multinomial coefficient for each element of the sample space
#'
#' @returns A numeric vector of probabilities
#'
#' @export
#'

calc_prob_null <- function(theta_cands, psi, psi0, SSpacearr, logC, II) {


  n <- sum(II)
  m <- ncol(SSpacearr)
  SSpacearr <- SSpacearr[II,]
  logC <- logC[II]

  res <- rep(NA, nrow(theta_cands))

  for(i in 1:nrow(theta_cands)) {

    thistheta <- theta_cands[i,]

      res[i] <- sum(exp((.colSums(t(SSpacearr) * log(thistheta), m = m, n = n) +
                           logC))) ## way faster
      #res[i] <- sum(exp((c(SSpacearr %*% log(thistheta)) + logC)[II]))

  }

  res[!is.na(res)]

}


#' Gradient of the multinomial likelihood sum
#'
#' @param theta_cands A matrix with samples in the rows and the parameters in
#'   the columns
#' @param psi The function of interest mapping parameters to the real line
#' @param psi0 The null boundary for testing psi <= psi0
#' @param minus1 Either plus or minus 1
#' @param SSpacearr A matrix with the sample space for the given size of the
#'   problem
#' @param II logical vector of sample space psi being more extreme than the
#'   observed psi
#'
#' @returns A matrix the same dimension as theta_cands
#'
#' @export
#' @examples
#' calc_prob_null_gradient(t(c(.28, .32, .4)),
#' \(s) s[, 1], .3, 1, matrix(c(2, 2, 1, 1, 2, 2, 0, 3, 2), ncol = 3),
#' rep(TRUE, 3))
#'
#' # numerically
#' testenv <- new.env()
#' testenv$SSpacearr <- matrix(c(2, 2, 1, 1, 2, 2, 0, 3, 2), ncol = 3)
#' testenv$thistheta <- c(.28, .32, .4)
#' numericDeriv(quote(sum(exp((.colSums(t(SSpacearr) * log(thistheta), m = 3, n = 3))))),
#'     theta = "thistheta", rho = testenv, central = TRUE)
#'
#'

calc_prob_null_gradient <- function(theta_cands, psi, psi0, SSpacearr, II) {

  m <- nrow(SSpacearr)
  n <- ncol(SSpacearr)
  SSpacearr <- SSpacearr[II,]

  res <- matrix(NA, nrow = nrow(theta_cands), ncol = ncol(theta_cands))

  for(i in 1:nrow(theta_cands)) {

    thistheta <- theta_cands[i,]

    res[i,] <- .colSums(t(t(SSpacearr) / thistheta) *
                         exp((.colSums(t(SSpacearr) * log(thistheta), m = n, n = sum(II)))), m = sum(II),
                        n = length(thistheta)) ## way faster


  }

  res


}

#' Calculate probability for given parameters
#'
#' Given a set of candidate parameter vectors, check if the null \eqn{\psi \leq
#' \psi_0} is satisfied, and if so, compute the probability for each element of
#' the sample space
#'
#' @param theta_cands A matrix with samples in the rows and the parameters in
#'   the columns
#' @param psi The function of interest mapping parameters to the real line
#' @param psi0 The null boundary for testing psi <= psi0
#' @param SSpacearr A matrix with the sample space for the given size of the
#'   problem
#' @param II logical vector of sample space psi being more extreme than the
#'   observed psi
#' @param logC log multinomial coefficient for each element of the sample space
#' @param psi_v Is psi vectorized by row?
#'
#' @returns A numeric vector of probabilities
#'
#' @export
#'

calc_prob_null2 <- function(theta_cands, psi, psi0, SSpacearr, logC, II) {

  SSpacearr <- SSpacearr[II,, drop = FALSE]
  logC <- logC[II]

  res <- calc_probs_rust(1.0*c(t(SSpacearr)), log(c(t(theta_cands))),
                         logC, d = ncol(SSpacearr), n = sum(II), nt = nrow(theta_cands))


  res

}



