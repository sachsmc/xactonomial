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

calc_prob_null <- function(theta_cands, psi, psi0, minus1, SSspacearr, logC, II) {


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

xactonomial <- function(psi, data, alpha = .05, psi_limits,
                        maxit = 50, chunksize = 500) {

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
                          lower = TRUE, target = alpha / 2, psi_limits = psi_limits,
                          SSpacearr = SSpacearr, logC = logC) {

    minus1 <- if(lower) 1 else -1
    II <- if(lower) psi_hat >= psi_obs else psi_hat <= psi_obs

    seqmaxes <- rep(NA, maxit)
    for(i in 1:maxit) {
      theta_cands <- do.call(cbind, lapply(d_k, \(i) get_theta_random(i, chunksize)))
      #theta_cands <- smart_sample_theta(psi, psi0 = psi0,
      #                                  psi_limits = psi_limits, d_k = d_k,
      #                                  chunksize = chunksize, lower = lower)
      these_probs <- calc_prob_null(theta_cands, psi, psi0, minus1,
                                    SSpacearr, logC, II)
      if(length(these_probs) == 0) next

      cand <- c(seqmaxes, these_probs)
      if(all(is.na(cand))) seqmaxes[i] <- 1e-12 else {
        seqmaxes[i] <- max(cand, na.rm = TRUE)
      }
      if(seqmaxes[i] > target + .001) break

    }
    if(all(is.na(seqmaxes))) return(1e-12) else max(seqmaxes, na.rm = TRUE)
  }


  flower <- function(x){
    pvalue_psi0(x, maxit = maxit, chunksize = chunksize,
          lower = TRUE, target = alpha / 2, psi_limits = psi_limits,
          SSpacearr = SSpacearr, logC = logC) - alpha / 2
  }

  fupper <- function(x) {
    pvalue_psi0(x, maxit = maxit, chunksize = chunksize,
                lower = FALSE, target = alpha / 2, psi_limits = psi_limits,
                SSpacearr = SSpacearr, logC = logC) - (alpha / 2)
  }



  lower_limit <- itp_root(flower, psi_limits[1], psi_limits[2],
                          fa = -alpha / 2, fb = 1 - alpha / 2, maxit = 10)
  upper_limit <- itp_root(fupper, psi_limits[1], psi_limits[2],
                          fa = 1-alpha / 2, fb = - alpha / 2, maxit = 10)


  list(estimate = psi_obs,
    conf.int = c(lower_limit, upper_limit),
    pvalue_function = pvalue_psi0)



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



itp_root <- function(f, a, b, k1 = .1, k2 = 2, n0 = 1,
                     eps = .005, maxit = 100, fa = NULL, fb = NULL,
                     verbose = FALSE) {

  if(is.null(fa) | is.null(fb)) {
    fa <- f(a)
    fb <- f(b)
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

    yitp <- f(xitp)
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

    if((b - a) < 2 * eps | j >= maxit) break

    for_rk <- for_rk * .5
    j <- j + 1

  }

  (a + b) / 2


}

