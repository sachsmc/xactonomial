#' Exact inference for a real-valued function of multinomial parameters
#'
#' We consider the k sample multinomial problem where we observe k vectors
#' (possibly of different lengths), each representing an independent sample from
#' a multinomial. For a given function psi which takes in the concatenated
#' vector of multinomial probabilities and outputs a real number, we are
#' interested in computing a p-value for a test of psi >= psi0, and constructing
#' a confidence interval for psi.
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
#' The computation is somewhat involved so it is best for small sample sizes
#' though there are some strategies for speeding things up that are described in Details.
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
#' @param psi.is.vectorized Logical. If TRUE, expect that psi can take a matrix as input, and return a vector of length the number of rows, computing the statistic for each row of the matrix. If possible, this will substantially speed up the computation. See examples.
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
#' xactonomial(psi_ba_v, data, psi_limits = c(0, 1), maxit = 5, chunksize = 20, psi.is.vectorized = TRUE)
#'

xactonomial <- function(data, psi, statistic = NULL, psi0 = NULL,
                        alternative = c("two.sided", "less", "greater"),
                        psi_is_vectorized = FALSE, psi_limits, theta_boundary_points = NULL, p_target = 1,
                        conf_int = TRUE, conf_level = .95, itp_maxit = 10,
                        maxit = 50, chunksize = 500,
                        theta_sampler = runif_dk_vects,
                        ga = TRUE, ga_gfactor = 1, ga_lrate = .01,
                        ga_restart_every = 10
                        ) {

  alternative <- match.arg(alternative)

  alpha <- 1 - conf_level

  k <- length(data)
  d_k <- sapply(data, length)
  n_k <- sapply(data, sum)

  sspace_size <- prod(sapply(data, \(dd) choose(length(dd) + sum(dd) - 1, length(dd) - 1)))

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


  if(is.null(statistic)) { ## compute psi with empirical proportions

    statistic <- function(df) {
      denom <- rep.int(n_k, d_k)
      if(is.matrix(df)) {
        pmat <- matrix(rep(denom, nrow(df)), nrow = nrow(df), ncol = length(denom), byrow = TRUE)
        if(psi_is_vectorized) {
          psi(df / pmat)
        } else {
          apply(df / pmat, 1, psi)
        }
      } else {
        if(psi_is_vectorized) {
          psi(t(df / denom))
        } else {
          psi(df / denom)
        }

      }

    }

  }

  psi_obs <- statistic(unlist(data))
  # psi_obs <- if(psi.is.vectorized) {
  #   psi(matrix(dprobs, nrow = 1))
  # } else {
  #   psi(tmpdat)
  # }

  psi_hat <- statistic(spacelist$Sspace)

  ## check if all elements of sample space above or below observed





  pvalues <- if(!is.null(psi0)) {

    if(!is.null(theta_boundary_points) & any(abs(psi0 - psi_limits) < 1e-8)) {

      psi0bnd <- if(all.equal(psi0, psi_limits[1], check.names = FALSE, check.attributes = FALSE)) {
        theta_boundary_points$lower
      } else {
        theta_boundary_points$upper
      }
      II.lower <- psi_hat >= psi_obs
      II.upper <- psi_hat <= psi_obs
     c(calc_prob_null(psi0bnd, spacelist$Sspace, spacelist$logC,  II.lower),
       calc_prob_null(psi0bnd, spacelist$Sspace, spacelist$logC,  II.upper))


    } else {

    pvalue_psi0(psi0 = psi0, psi = psi, psi_hat = psi_hat, psi_obs = psi_obs, alternative = alternative,
                maxit = maxit, chunksize = chunksize,
                            p_target = p_target, SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                            d_k = d_k, psi_is_vectorized = psi_is_vectorized,
                            theta_sampler = theta_sampler,
                            ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                ga_restart_every = ga_restart_every)

    }

    } else NA


  confint <- if(conf_int) {
  flower <- function(x, psi, psi_hat, psi_obs, maxit, chunksize, p_target, SSpacearr, logC, d_k, psi_is_vectorized,
                     theta_sampler, ga, ga_gfactor, ga_lrate, ga_restart_every){
    pvalue_psi0(x, psi, psi_hat, psi_obs, maxit = maxit, chunksize = chunksize,
                p_target = p_target, SSpacearr = SSpacearr, logC = logC,
                d_k = d_k, psi_is_vectorized = psi_is_vectorized,
                theta_sampler = theta_sampler,
                ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                ga_restart_every = ga_restart_every)[1] - alpha / 2}

  if(!is.null(theta_boundary_points$lower)) {

    II.lower <- psi_hat >= psi_obs

    flower.boundary <- max(calc_prob_null2(theta_boundary_points$lower, spacelist$Sspace, spacelist$logC,  II.lower)) -
      alpha / 2

  } else {
    flower.boundary <- -alpha / 2
  }


  fupper <- function(x, psi, psi_hat, psi_obs, maxit, chunksize, p_target, SSpacearr, logC, d_k, psi_is_vectorized,
                     theta_sampler, ga, ga_gfactor, ga_lrate, ga_restart_every){
    pvalue_psi0(x, psi, psi_hat, psi_obs, maxit = maxit, chunksize = chunksize,
                p_target = p_target, SSpacearr = SSpacearr, logC = logC,
                d_k = d_k, psi_is_vectorized = psi_is_vectorized,
                theta_sampler = theta_sampler,
                ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                ga_restart_every = ga_restart_every)[2] - alpha / 2
  }

  if(!is.null(theta_boundary_points$upper)) {

    II.upper <- psi_hat <= psi_obs
    fupper.boundary <- max(calc_prob_null2(theta_boundary_points$upper, spacelist$Sspace, spacelist$logC,  II.upper)) -
      alpha / 2

  } else {
    fupper.boundary <- - alpha / 2
  }


  if(isTRUE(all.equal(psi_obs, min(psi_hat))) | flower.boundary > 0) {
    lower_limit <- psi_limits[1]
  } else {
    lower_limit <- itp_root(flower, psi_limits[1], psi_limits[2],
                            fa = flower.boundary, fb = 1 - alpha / 2, maxiter = itp_maxit,
                            psi = psi, psi_hat = psi_hat, psi_obs = psi_obs,
                            maxit = maxit, chunksize = chunksize,
                            p_target = alpha / 2 + .01,
                            SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                            d_k = d_k, psi_is_vectorized = psi_is_vectorized,
                            theta_sampler = theta_sampler,
                            ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                            ga_restart_every = ga_restart_every)

    if((attr(lower_limit, "iter") %||% itp_maxit) == itp_maxit){
      lower_limit <- psi_limits[1]
      warning("No root found for lower confidence limit, using lower boundary.")
    }

  }


  if(isTRUE(all.equal(psi_obs, max(psi_hat))) | fupper.boundary > 0) {
    upper_limit <- psi_limits[2]
  } else {
    upper_limit <- itp_root(fupper, psi_limits[1], psi_limits[2],
                            fa = 1-alpha / 2, fb = fupper.boundary, maxiter = itp_maxit,
                            psi = psi, psi_hat = psi_hat, psi_obs = psi_obs,
                            maxit = maxit, chunksize = chunksize,
                            p_target = alpha / 2 + 0.01,
                            SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                            d_k = d_k, psi_is_vectorized = psi_is_vectorized,
                            theta_sampler = theta_sampler,
                            ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                            ga_restart_every = ga_restart_every)

    if((attr(upper_limit, "iter") %||%  itp_maxit) == itp_maxit){
      upper_limit <- psi_limits[2]
      warning("No root found for upper confidence limit, using upper boundary.")
    }

  }


  c(lower_limit, upper_limit)
  } else c(NA, NA)

  p.sequence <- attr(pvalues, "p.sequence")

  # list(estimate = psi_obs,
  #      conf.int = confint,
  #      p.value = c(pvalues, 2 * min(pvalues)),
  #      p.sequence = p.sequence
  #      )

  attr(confint, "conf.level") = 1 - alpha

  res <- list(
    estimate = psi_obs,
    p.value = switch(alternative, "greater" = pvalues[1], "less" = pvalues[2], "two.sided" = 2 * min(pvalues)),
    conf.int = confint,
    null.value = c(psi0 = psi0),
    alternative = alternative,
    method = "Monte-Carlo exact multinomial test",
    data.name = deparse1(substitute(data)),
    p.sequence = p.sequence
  )
  class(res) <- "htest"
  res

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
#' @param theta_sampler Sampler
#' @param gradient_ascent Do ga?
#' @param gamma Parameter to determine concentration of dirichlet, only relevant if gradient_ascent = TRUE
#' @returns A p-value
#'

pvalue_psi0 <- function(psi0, psi, psi_hat, psi_obs, alternative = "two.sided",
                        maxit, chunksize,
                        p_target,
                        SSpacearr, logC, d_k, psi_is_vectorized = FALSE,
                        theta_sampler = runif_dk_vects,
                        ga = FALSE, ga_gfactor = 1, ga_lrate = .01,
                        ga_restart_every = 10
                        ) {

  II.lower <- psi_hat >= psi_obs
  II.upper <- psi_hat <= psi_obs

  p.null <- rep(1 / (maxit * chunksize), maxit)
  p.alt <- rep(1 / (maxit * chunksize), maxit)

  theta_cands <- theta_sampler(d_k, chunksize)
  null_continue <- alt_continue <- TRUE
  if(alternative == "greater") alt_continue <- FALSE
  if(alternative == "lower") null_continue <- FALSE
  null_stop <- alt_stop <- maxit
  never_null <- null_continue
  never_alt <- alt_continue
  for(i in 1:maxit) {

    if(isFALSE(null_continue) & isFALSE(alt_continue)) break

    this_theta <- theta_cands
    psi_theta <- if(psi_is_vectorized) psi(theta_cands) else apply(theta_cands, MAR = 1, psi)
    null_indicator <- psi_theta <= psi0

    if(i %% ga_restart_every == 0) {
      theta_cands <- theta_sampler(d_k, chunksize)
      p.null[i] <- p.null[i-1]
      p.alt[i] <- p.alt[i-1]
      next
    }

    if(sum(null_indicator) > 0 & isTRUE(null_continue)) {

      never_null <- FALSE
      theta_null <- this_theta[null_indicator, , drop = FALSE]
      probs_null <- calc_prob_null2(theta_null,
                                    SSpacearr, logC, II.lower)

      p.null[i] <- max(c(p.null, probs_null), na.rm = TRUE)
      if(p.null[i] > p_target) {
        null_continue <- FALSE
        null_stop <- i
      }

      if(isTRUE(null_continue)) {
      if(ga) {
        grad_this <- calc_prob_null_gradient(theta_null[which.max(probs_null), , drop = FALSE],
                                             SSpacearr, II.lower)

        theta_cands_n <- theta_null[which.max(probs_null), , drop = FALSE] +
          ga_lrate * grad_this[1,]
        theta_cands_n <- theta_cands_n / sum(theta_cands_n)
        gammat <- if(ga_gfactor == "adapt") 1 / max(.001, min(theta_cands_n)) else ga_gfactor
        theta_cands_n <- rdirich_dk_vects(chunksize, list(gammat * theta_cands_n))
      } else {
        theta_cands_n <- theta_sampler(d_k, chunksize)
      }
      }

    } else {
      if(isTRUE(null_continue)) {
      theta_cands_n <- theta_sampler(d_k, chunksize)
      p.null[i] <- if(i > 1) p.null[i-1] else 1 / (maxit * chunksize)
      }

    }

    if(sum(!null_indicator) > 0 & isTRUE(alt_continue)) {
      never_alt <- FALSE
      theta_alt <- this_theta[!null_indicator, , drop = FALSE]
      probs_alt <- calc_prob_null2(theta_alt,
                                 SSpacearr, logC, II.upper)
    p.alt[i] <- max(c(p.alt, probs_alt), na.rm = TRUE)
    if(p.alt[i] > p_target) {
      alt_continue <- FALSE
      alt_stop <- i
    }

    if(isTRUE(alt_continue)){
    ## check the gradient at the largest value
    if(ga) {

      grad_this_alt <- calc_prob_null_gradient(theta_alt[which.max(probs_alt), , drop = FALSE],
                                            SSpacearr, II.upper)
      theta_cands_a <- theta_alt[which.max(probs_alt), , drop = FALSE] + ga_lrate * grad_this_alt[1,]
      theta_cands_a <- theta_cands_a / sum(theta_cands_a)
      gammat <- if(ga_gfactor == "adapt") 1 / max(.001, min(theta_cands_a)) else ga_gfactor
      theta_cands_a <- rdirich_dk_vects(chunksize, list(gammat * theta_cands_a))

    } else {
      theta_cands_a <- theta_sampler(d_k, chunksize)
    }
    }

    } else {
      if(isTRUE(alt_continue)) {
      theta_cands_a <- theta_sampler(d_k, chunksize)
      p.alt[i] <- if(i > 1) p.alt[i-1] else 1 / (maxit * chunksize)
      }
    }

    theta_cands <- if(null_continue & alt_continue) {
      rbind(theta_cands_n, theta_cands_a)
    } else if(null_continue){
      theta_cands_n
    } else if(alt_continue) {
      theta_cands_a
    }


  }
  if(isTRUE(never_null) | isTRUE(never_alt)) {
    warning("Never found a value from the null space, check results carefully and consider increasing iterations!")
  }

  res <- c(null = max(p.null, na.rm = TRUE), alt = max(p.alt, na.rm = TRUE))
  attr(res, "p.sequence") <- list(p.null = p.null[1:null_stop], p.alt = p.alt[1:alt_stop])
  res
}
