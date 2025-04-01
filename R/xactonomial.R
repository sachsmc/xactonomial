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
#' \eqn{H_0: \psi(\boldsymbol{\theta}) \neq \psi_0} for the two sided case,
#' \eqn{H_0: \psi(\boldsymbol{\theta}) \leq \psi_0} for the case \code{alternative = "greater"}, and
#' \eqn{H_0: \psi(\boldsymbol{\theta}) \geq \psi_0} for the case \code{alternative = "less"}.
#' We make no assumptions and do not rely on large sample approximations.
#' The computation is somewhat involved so it is best for small sample sizes.
#'
#' @param psi Function that takes in a vector of parameters and outputs a real
#'   valued number
#' @param data A list with k elements representing the vectors of counts of a
#'   k-sample multinomial
#' @param psi0 The null hypothesis value for the parameter being tested. A p value for a test of psi <= psi0 is computed. If NULL only a confidence interval is computed.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param alpha A 1 - alpha percent confidence interval will be computed
#' @param psi_limits A vector of length 2 giving the lower and upper limits of
#'   the range of \eqn{\psi(\theta)}
#' @param maxit Maximum number of iterations of the stochastic procedure
#' @param chunksize The number of samples taken from the parameter space at each
#'   iteration
#' @param conf.int Logical. If FALSE, no confidence interval is calculated, only the p-value.
#' @param psi_is_vectorized Logical. If TRUE, expect that psi can take a matrix as input, and return a vector of length the number of rows, computing the statistic for each row of the matrix. If possible, this will substantially speed up the computation. See examples.
#'
#' @returns A list with 3 elements: the estimate, the 1 - alpha percent
#'   confidence interval, and p-value
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
#' psi_ba_v <- function(theta) {
#' theta1 <- theta[,1:4, drop = FALSE]
#' theta2 <- theta[,5:8, drop = FALSE]
#' rowSums(sqrt(theta1 * theta2))
#' }
#' data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
#' xactonomial(psi_ba_v, data, psi_limits = c(0, 1), maxit = 5, chunksize = 20, psi_is_vectorized = TRUE)
#'

xactonomial <- function(psi, data, psi0 = NULL, alternative = c("two.sided", "less", "greater"),
                        alpha = .05, psi_limits,
                        maxit = 50, chunksize = 500, conf.int = TRUE,
                        psi_is_vectorized = FALSE
                        ) {

  alt <- match.arg(alternative)

  k <- length(data)
  d_k <- sapply(data, length)

  tmpdat <- lapply(data, \(x) x / sum(x)) |> unlist()
  psi_obs <- if(psi_is_vectorized) {
    psi(matrix(tmpdat, nrow = 1))
  } else {
    psi(tmpdat)
  }

  SSpace <- lapply(data, \(x) sspace_multinom(length(x), sum(x)))

  if(k == 1) {

    newX <- SSpace[[1]]
    sumX <- sum(newX[1,])
    logC <- lfactorial(sumX) - rowSums(apply(newX, 2, lfactorial))

    spacelist <- list(Sspace = newX, Sprobs = newX / sumX, logC = logC)

  } else if(k == 2) {

    spacelist <- combinate(SSpace[[1]], SSpace[[2]])

  } else {

    spacelist <- combinate(SSpace[[1]], SSpace[[2]])
    for(i in 3:k) {

      spacelist <- combinate2(spacelist, SSpace[[i]])

    }

  }


  psi_hat <- if(psi_is_vectorized) psi(spacelist$Sprobs) else apply(spacelist$Sprobs, 1, psi)


  pvalue <- if(!is.null(psi0)) {

    if(alt == "greater") {
      pvalue_psi0(psi0 = psi0, psi = psi, psi_hat = psi_hat,
                           psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                           lower = TRUE, target = alpha / 2,
                           SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                  psi_v = psi_is_vectorized)
    } else if(alt == "less") {

      pvalue_psi0(psi0 = psi0, psi = psi, psi_hat = psi_hat,
                  psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                  lower = FALSE, target = alpha / 2,
                  SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                  psi_v = psi_is_vectorized)

    } else if(alt == "two.sided") {

      pl <- pvalue_psi0(psi0 = psi0, psi = psi, psi_hat = psi_hat,
                        psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                        lower = TRUE, target = alpha / 2,
                        SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                        psi_v = psi_is_vectorized)
      pu <- pvalue_psi0(psi0 = psi0, psi = psi, psi_hat = psi_hat,
                        psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                        lower = FALSE, target = alpha / 2,
                        SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                        psi_v = psi_is_vectorized)

      2 * min(pl, pu)

    }

  } else NA


  confint <- if(conf.int) {
  flower <- function(x, psi, psi_hat, psi_obs, maxit, chunksize, target, SSpacearr, logC, d_k, psi_v){
    pvalue_psi0(psi0 = x, psi = psi, psi_hat = psi_hat, psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                lower = TRUE, target = alpha / 2,
                SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                psi_v = psi_v) - alpha / 2
  }

  fupper <- function(x, psi, psi_hat, psi_obs, maxit, chunksize, target, SSpacearr, logC, d_k, psi_v) {
    pvalue_psi0(psi0 = x, psi = psi, psi_hat = psi_hat, psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                lower = FALSE, target = alpha / 2,
                SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                psi_v = psi_v) - (alpha / 2)
  }


  if(isTRUE(all.equal(psi_obs, psi_limits[1]))) {
    lower_limit <- psi_obs
  } else {
    lower_limit <- itp_root(flower, psi_limits[1], psi_limits[2],
                            fa = -alpha / 2, fb = 1 - alpha / 2, maxiter = 10,
                            psi = psi, psi_hat = psi_hat, psi_obs = psi_obs,
                            maxit = maxit, chunksize = chunksize,
                            target = alpha / 2,
                            SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                            psi_v = psi_is_vectorized)

  }


  if(isTRUE(all.equal(psi_obs, psi_limits[2]))) {
    upper_limit <- psi_obs
  } else {
    upper_limit <- itp_root(fupper, psi_limits[1], psi_limits[2],
                            fa = 1-alpha / 2, fb = - alpha / 2, maxiter = 10,
                            psi = psi, psi_hat = psi_hat, psi_obs = psi_obs,
                            maxit = maxit, chunksize = chunksize,
                            target = alpha / 2,
                            SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                            psi_v = psi_is_vectorized)

  }


  c(lower_limit, upper_limit)
  } else c(NA, NA)

  list(estimate = psi_obs,
       conf.int = confint,
       p.value = pvalue
       )



}


#' Compute a p value for the test of psi <= psi0 (lower = TRUE) or psi >= psi0 (lower = FALSE)
#'
#' @param psi0 The null value
#' @param psi The function of interest
#' @param psi_hat The vector of psi values at each element of the sample space
#' @param psi_obs The observed estimate
#' @param maxit Maximum iterations
#' @param chunksize Chunk size
#' @param lower Do a one sided test of the null that it is less than psi0, otherwise greater.
#' @param target Stop the algorithm if p >= target (for speed)
#' @param SSpacearr The sample space array
#' @param logC The log multinomial coefficient
#' @param d_k The vector of dimensions
#' @param psi_v Is psi vectorized by row?
#' @returns A p-value
#'

pvalue_psi0 <- function(psi0, psi, psi_hat, psi_obs, maxit, chunksize,
                        lower = TRUE, target,
                        SSpacearr, logC, d_k, psi_v = FALSE,
                        sample_theta = runif_dk_vects
                        ) {

  minus1 <- if(lower) 1 else -1
  II <- if(lower) psi_hat >= psi_obs else psi_hat <= psi_obs

  seqmaxes <- rep(1 / (maxit * chunksize), maxit)
  for(i in 1:maxit) {
    theta_cands <- sample_theta(d_k, chunksize) #do.call("cbind", lapply(d_k, \(i) get_theta_random(i, chunksize)))

    these_probs <- calc_prob_null2(theta_cands, psi, psi0, minus1,
                                   SSpacearr, logC, II, psi_v = psi_v)

    if(length(these_probs) == 0) next

    cand <- c(seqmaxes, these_probs)
    if(all(is.na(cand))) seqmaxes[i] <- 1/(maxit * chunksize) else {
      seqmaxes[i] <- max(cand, na.rm = TRUE)
    }
    if(seqmaxes[i] > target + .001) break

  }
  max(seqmaxes, na.rm = TRUE)
}
