
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


  flower <- function(x, maxit, chunksize, target, psi_limits, SSpacearr, logC){
    pvalue_psi0(x, maxit = maxit, chunksize = chunksize,
                lower = TRUE, target = alpha / 2, psi_limits = psi_limits,
                SSpacearr = SSpacearr, logC = logC) - alpha / 2
  }

  fupper <- function(x, maxit, chunksize, target, psi_limits, SSpacearr, logC) {
    pvalue_psi0(x, maxit = maxit, chunksize = chunksize,
                lower = FALSE, target = alpha / 2, psi_limits = psi_limits,
                SSpacearr = SSpacearr, logC = logC) - (alpha / 2)
  }



  lower_limit <- itp_root(flower, psi_limits[1], psi_limits[2],
                          fa = -alpha / 2, fb = 1 - alpha / 2, maxiter = 10,
                          maxit = maxit, chunksize = chunksize,
                          target = alpha / 2, psi_limits = psi_limits,
                          SSpacearr = SSpacearr, logC = logC)
  upper_limit <- itp_root(fupper, psi_limits[1], psi_limits[2],
                          fa = 1-alpha / 2, fb = - alpha / 2, maxiter = 10,
                          maxit = maxit, chunksize = chunksize,
                          target = alpha / 2, psi_limits = psi_limits,
                          SSpacearr = SSpacearr, logC = logC)


  list(estimate = psi_obs,
       conf.int = c(lower_limit, upper_limit),
       pvalue_function = pvalue_psi0)



}
