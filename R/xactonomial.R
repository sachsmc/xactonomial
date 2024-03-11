#' Exact inference for a real-valued function of multinomial parameters
#'
#' We consider the k sample multinomial problem where we observe k vectors
#' (possibly of different lengths), each representing an independent sample from
#' a multinomial. For a given function psi which takes in the concatenated
#' vector of multinomial probabilities and outputs a real number, we are
#' interested in constructing a confidence interval for psi.
#'
#' Let \eqn{T_j} be distributed
#' \eqn{\mbox{Multinomial}_{d_j}(\boldsymbol{\theta}_j, n_j)} for \eqn{j = 1,
#' \ldots, k} and denote \eqn{\boldsymbol{T} = (T_1, \ldots, T_k)} and
#' \eqn{\boldsymbol{\theta} = (\theta_1, \ldots, \theta_k)}. The subscript
#' \eqn{d_j} denotes the dimension of the multinomial. Suppose one is interested
#' in the parameter \eqn{\psi(\boldsymbol{\theta}) \in \Theta \subseteq
#' \mathbb{R}}. Given a sample of size \eqn{n} from \eqn{\boldsymbol{T}}, one
#' can estimate \eqn{\boldsymbol{\theta}} with the sample proportions as
#' \eqn{\hat{\boldsymbol{\theta}}} and hence
#' \eqn{\pi(\hat{\boldsymbol{\theta}})}. This function constructs a \eqn{1 -
#' \alpha} percent confidence interval for \eqn{\psi(\boldsymbol{\theta})} and
#' provides a function to calculate a p value for a test of the null hypothesis
#' \eqn{H_0: \psi(\boldsymbol{\theta}) \leq \psi_0}. We make no assumptions and
#' do not rely on large sample approximations. The computation in somewhat
#' involved so it is best for small sample sizes.
#'
#' @param psi Function that takes in a vector of parameters and outputs a real
#'   valued number
#' @param data A list with k elements representing the vectors of counts of a
#'   k-sample multinomial
#' @param alpha A 1 - alpha percent confidence interval will be computed
#' @param psi_limits A vector of length 2 giving the lower and upper limits of
#'   the range of \eqn{\psi(\theta)}
#' @param maxit Maximum number of iterations of the stochastic procedure
#' @param chunksize The number of samples taken from the parameter space at each
#'   iteration
#'
#' @returns A list with 3 elements: the estimate, the 1 - alpha percent
#'   confidence interval, and a function for calculation of p-values
#' @export
#' @examples
#' psi_ba <- function(theta) {
#'   theta1 <- theta[1:4]
#'   theta2 <- theta[5:8]
#'   sum(sqrt(theta1 * theta2))
#'   }
#' data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
#' xactonomial(psi_ba, data, psi_limits = c(0, 1), maxit = 5, chunksize = 20)
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
      theta_cands <- do.call("cbind", lapply(d_k, \(i) get_theta_random(i, chunksize)))

      these_probs <- calc_prob_null2(theta_cands, psi, psi0, minus1,
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
