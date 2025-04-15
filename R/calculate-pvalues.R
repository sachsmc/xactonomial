#' @inherit calc_prob_null_fast
#' @export
#'

calc_prob_null <- function(theta_cands, SSpacearr, logC, II) {


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
#' matrix(c(2, 2, 1, 1, 2, 2, 0, 3, 2), ncol = 3),
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

calc_prob_null_gradient <- function(theta_cands, SSpacearr, II) {

  m <- nrow(SSpacearr)
  n <- ncol(SSpacearr)
  SSpacearr <- SSpacearr[II,]

  theta_cands[theta_cands < 1e-250] <- 1e-250
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
#' Given a set of candidate parameter vectors, the enumerated sample space, and a
#' logical vector with the same number of elements of the sample space, compute
#' the probability for each element of the sample space and take the sum.
#'
#' @param theta_cands A matrix with samples in the rows and the parameters in
#'   the columns
#' @param SSpacearr A matrix with the sample space for the given size of the
#'   problem
#' @param II logical vector of sample space psi being more extreme than the
#'   observed psi
#' @param logC log multinomial coefficient for each element of the sample space
#'
#' @returns A numeric vector of probabilities
#'
#' @export
#' @examples
#' sspace_3_5 <- matrix(sspace_multinom(3, 5), ncol = 3, byrow = TRUE)
#' theta_cands <- matrix(sample_unit_simplexn(3, 10), ncol = 3,byrow = TRUE)
#' calc_prob_null_fast(theta_cands, sspace_3_5,
#' apply(sspace_3_5, 1, log_multinom_coef, sumx = 5), II = 1:21 > 12)
#' # same as below but faster
#' calc_prob_null(theta_cands, sspace_3_5,
#' apply(sspace_3_5, 1, log_multinom_coef, sumx = 5), II = 1:21 > 12)
#'

calc_prob_null_fast <- function(theta_cands, SSpacearr, logC, II) {

  SSpacearr <- SSpacearr[II,, drop = FALSE]
  logC <- logC[II]

  theta_cands[theta_cands < 1e-250] <- 1e-250
  res <- calc_multinom_probs(1.0*c(t(SSpacearr)), log(c(t(theta_cands))),
                         logC, d = ncol(SSpacearr), n = sum(II), nt = nrow(theta_cands))


  res

}



