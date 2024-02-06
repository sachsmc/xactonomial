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


  checkpsi <- minus1 * apply(theta_cands, MAR = 1, psi) <= minus1 * psi0

  if(sum(checkpsi) == 0) return(NA)
  theta_cands <- theta_cands[checkpsi, , drop = FALSE]


  res <- rep(NA, nrow(theta_cands))

  for(i in 1:nrow(theta_cands)) {

    thistheta <- theta_cands[i,]
    thispsi <- psi(thistheta)
    if(minus1 * thispsi <= minus1 * psi0) {

      res[i] <- sum(exp((colSums(t(SSpacearr) * log(thistheta)) + logC)[II])) ## way faster
      #res[i] <- sum(exp((c(SSpacearr %*% log(thistheta)) + logC)[II]))

    }

  }

  res[!is.na(res)]

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


#' Do the damn thing
#'
#' It is assumed that the parameters given to psi are in the same order as given in d_k
#'
#' @param psi Function that takes in a vector of parameters and outputs a real valued number
#' @param data A list with k elements representing the vectors of counts of a k-sample multinomial
#' @param alpha A 1 - alpha percent confidence interval will be computed
#'

xactonomial <- function(psi, data, alpha = .05,
                        maxit = 500, chunksize = 500) {

  k <- length(data)
  d_k <- sapply(data, length)
  SSpace <- lapply(data, \(x) sspace_multinom(length(x), sum(x)))

  bigdex <- expand_index(sapply(SSpace, nrow))
  psi_hat <- logC <- rep(NA, nrow(bigdex))
  psi_obs <- do.call(psi, lapply(data, \(x) x / sum(x)) |> unlist() |> list())

  SSpacearr <- array(dim = c(nrow(bigdex), sum(d_k)))

  for(i in 1:nrow(bigdex)) {

    thisS <- lapply(1:length(bigdex[i,]), \(j){
      Sj <- SSpace[[j]][bigdex[i, j],]
      Sj
    })
    SSpacearr[i,] <- unlist(thisS)
    psi_hat[i] <- do.call(psi, lapply(thisS, \(x) x / sum(x)) |> unlist() |> list())
    logC[i] <- sum(sapply(thisS, \(x) log_multinom_coef(x, sum(x))))

  }

  pvalue_psi0 <- function(psi0, maxit = maxit, chunksize = chunksize,
                          lower = TRUE, target = alpha / 2,
                          SSpacearr = SSpacearr, logC = logC) {

    minus1 <- if(lower) 1 else -1
    II <- if(lower) psi_hat >= psi_obs else psi_hat <= psi_obs

    seqmaxes <- rep(NA, maxit)
    for(i in 1:maxit) {
      theta_cands <- do.call(cbind, lapply(d_k, \(i) get_theta_random(i, chunksize)))
      these_probs <- calc_prob_null(theta_cands, psi, psi0, minus1, SSpacearr, II, logC)
      if(length(these_probs) == 0) next

      seqmaxes[i] <- max(c(seqmaxes, these_probs), na.rm = TRUE)
      if(seqmaxes[i] > target + .001) break

    }
    max(seqmaxes, na.rm = TRUE)
  }


  lower_limit <- uniroot(\(x) pvalue_psi0(x, maxit = maxit, chunksize = chunksize, lower = TRUE, target = alpha/2,
                                          SSpacearr = SSpacearr, logC = logC) - alpha / 2,
                         f.lower = -alpha/2, f.upper = 1 - alpha/2,
                         interval = c(-100, 100), tol = .001)
  upper_limit <- uniroot(\(x) pvalue_psi0(x, maxit = maxit, chunksize = chunksize, lower = FALSE, target = alpha/2,
                                          SSpacearr = SSpacearr, logC = logC) - alpha / 2,
                         f.lower = 1 - alpha/2, f.upper = -alpha/2,
                         interval = c(-100, 100), tol = .001)


  list(estimate = psi_obs,
    conf.int = c(lower_limit$root, upper_limit$root),
    pvalue_function = pvalue_psi0)



}





