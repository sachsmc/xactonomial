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

calc_prob_null <- function(theta_cands, psi, psi0, minus1, SSpacearr, logC, II) {


  checkpsi <- if(minus1 == 1) {
    apply(theta_cands, MAR = 1, psi) <= psi0
  } else {
    apply(theta_cands, MAR = 1, psi) >= psi0
  }

  if(sum(checkpsi) == 0) return(NA)
  theta_cands <- theta_cands[checkpsi, , drop = FALSE]
  m <- nrow(SSpacearr)
  n <- ncol(SSpacearr)
  SSpacearr <- SSpacearr[II,]
  logC <- logC[II]

  res <- rep(NA, nrow(theta_cands))

  for(i in 1:nrow(theta_cands)) {

    thistheta <- theta_cands[i,]

      res[i] <- sum(exp((.colSums(t(SSpacearr) * log(thistheta), m = n, n = sum(II)) +
                           logC))) ## way faster
      #res[i] <- sum(exp((c(SSpacearr %*% log(thistheta)) + logC)[II]))

  }

  res[!is.na(res)]

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

calc_prob_null2 <- function(theta_cands, psi, psi0, minus1, SSpacearr, logC, II) {


  checkpsi <- if(minus1 == 1) {
    apply(theta_cands, MAR = 1, psi) <= psi0
  } else {
    apply(theta_cands, MAR = 1, psi) >= psi0
  }

  if(sum(checkpsi) == 0) return(NA)
  theta_cands <- theta_cands[checkpsi, , drop = FALSE]
  m <- nrow(SSpacearr)
  n <- ncol(SSpacearr)
  SSpacearr <- SSpacearr[II,]
  logC <- logC[II]


  res <- calc_probs_rust(1.0*c(t(SSpacearr)), log(c(t(theta_cands))),
                         logC, d = ncol(SSpacearr), n = sum(II), nt = nrow(theta_cands))


  res

}



