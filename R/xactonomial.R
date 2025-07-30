#' Improved inference for a real-valued function of multinomial parameters
#'
#' We consider the k sample multinomial problem where we observe k vectors
#' (possibly of different lengths), each representing an independent sample from
#' a multinomial. For a given function \eqn{\tau(\theta)} which takes in the concatenated
#' vector of multinomial probabilities \eqn{\theta} and outputs a real number, we are
#' interested in computing a p-value for a test of \eqn{\tau(\theta) = \psi \geq \psi_0}, and constructing
#' a confidence interval for \eqn{\psi}.
#'
#' Let \eqn{T_j} be distributed
#' \eqn{\mbox{Multinomial}_{d_j}(\boldsymbol{\theta}_j, n_j)} for \eqn{j = 1,
#' \ldots, k} and denote \eqn{\boldsymbol{T} = (T_1, \ldots, T_k)} and
#' \eqn{\boldsymbol{\theta} = (\theta_1, \ldots, \theta_k)}. The subscript
#' \eqn{d_j} denotes the dimension of the multinomial. Suppose one is interested
#' in the parameter \eqn{\psi = \tau(\boldsymbol{\theta}) \in \Psi \subseteq
#' \mathbb{R}}. Given a sample of size \eqn{n} from \eqn{\boldsymbol{T}}, say
#' \eqn{\boldsymbol{X} = (X_1, \ldots, X_k)}, which is a vector of counts obtained
#' by concatenating the k independent count vectors, let \eqn{G(\boldsymbol{X})}
#' denote a real-valued statistic that defines the ordering of the sample space.
#' The default choice of the statistic is to estimate \eqn{\boldsymbol{\theta}}
#' with the sample proportions and plug them into \eqn{\tau(\boldsymbol{\theta})}.
#' This function calculates a p value for a test of the null hypothesis
#' \eqn{H_0: \psi \neq \psi_0} for the two sided case,
#' \eqn{H_0: \psi \leq \psi_0} for the case \code{alternative = "greater"}, and
#' \eqn{H_0: \psi \geq \psi_0} for the case \code{alternative = "less"}.
#' We make no assumptions and do not rely on large sample approximations.
#' It also optionally constructs a \eqn{1 - \alpha} percent confidence interval for \eqn{\psi}.
#' The computation is somewhat involved so it is best for small sample sizes. The
#' calculation is done by sampling a large number of points from the null parameter space \eqn{\Theta_0},
#' then computing multinomial probabilities under those values for the range of the sample space
#' where the statistic is as or more extreme than the observed statistic given data. It
#' is basically the definition of a p-value implemented with Monte Carlo methods. Some
#' options for speeding up the calculation are available.
#'
#' @section Specifying the function \eqn{\tau(\cdot)}:
#' This function is the \code{f_param} argument and should be a function that either: 1)
#' takes a vector of length sum(d_j) (the total number of bins) and outputs a single number, or 2) takes a
#' matrix with number of columns equal to sum(d_j), and arbitrary number of rows and
#' outputs a vector with length equal to the number of rows. In other words, psi can be
#' not vectorized or it can be vectorized by row. Writing it so that it is vectorized
#' can speed up the calculation. See examples.
#'
#' @section Boundary issues:
#' It is required to provide \code{psi_limits}, a vector of length 2 giving the
#' smallest and largest possible values that the function psi can take, e.g., \code{c(0, 1)}.
#' If the null hypothesis value psi0 is at one of the limits, it is often the case
#' that sampling from the null parameter space is impossible because it is a set of
#' measure 0. While it may have measure 0, it is not empty, and will contain a finite
#' set of points. Thus you should provide the argument \code{theta_null_points} which is
#' a matrix where the rows contain the finite set (sometimes 1) of points
#' \eqn{\theta} such that \eqn{\tau(\theta) = \psi_0}. There is also an argument called
#' \code{p_value_limits} that can be used to improve performance of confidence intervals
#' around the boundary. This should be a vector of length 2 with the p-value for a test
#' of psi_0 <= psi_limits[1] and the p-value for a test of psi_0 >= psi_limits[2].
#' See examples.
#'
#' @section Optimization options:
#' For p-value calculation, you can provide a parameter p_target, so that the sampling
#' algorithm terminates when a p-value is found that exceeds p_target. The algorithm
#' begins by sampling uniformly from the unit simplices defining the parameter space, but
#' alternatives can be specified in \code{theta_sampler}. By default
#' gradient ascent (\code{ga = TRUE}) is performed during the p-value maximization
#' procedure, and \code{ga_gfactor} and \code{ga_lrate} control options for the gradient
#' ascent. At each iteration, the gradient of the multinomial probability at the current maximum
#' theta is computed, and a step is taken to \code{theta + lrate * gradient}. Then
#' for the next iteration, a set of \code{chunksize} samples are drawn from a Dirichlet distribution
#' with parameter \code{ga_gfactor * (theta + ga_lrate * gradient)}. If \code{ga_gfactor = "adapt"} then
#' it is set to \code{1 / max(theta)} at each iteration. The ITP algorithm \link{itp_root} is used
#' to find roots of the p-value function as a function of the psi0 value to get confidence intervals.
#' The maximum number of iterations and epsilon can be controlled via \code{itp_maxit, itp_eps}.
#'
#'
#'
#' @param data A list with k elements representing the vectors of counts of a
#'   k-sample multinomial
#' @param f_param Function that takes in parameters and outputs psi, a real
#'   valued number for each parameter. Can be vectorized rowwise for a matrix or not.
#' @param statistic Function that takes in a matrix with data vectors in the rows, and outputs a vector with the number of rows in the matrix. If NULL, will be inferred from \code{f_param} by plugging in the empirical proportions.
#' @param psi0 The null hypothesis value for the parameter being tested.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param psi_limits A vector of length 2 giving the lower and upper limits of
#'   the range of \eqn{\tau(\theta)}
#' @param theta_null_points An optional matrix where each row is a theta value that gives \code{f_param(theta) = psi0}. If this is supplied and psi0 = one of \code{psi_limits}, then a truly exact p-value will be calculated.
#' @param p_target If a p-value is found that is greater than p_target, terminate the algorithm early.
#' @param conf_int If TRUE, calculates a confidence interval by inverting the p-value function
#' @param conf_level A number between 0 and 1, the confidence level.
#' @param itp_maxit Maximum iterations to use in the ITP algorithm. Only relevant if conf_int = TRUE.
#' @param itp_eps Epsilon value to use for the ITP algorithm. Only relevant if conf_int = TRUE.
#' @param p_value_limits A vector of length 2 giving lower bounds on the one-sided p-values corresponding to psi0 = psi_limits, with alternative = "less" and "greater", respectively. Only relevant if conf_int = TRUE. See examples.
#' @param maxit Maximum number of iterations of the Monte Carlo procedure
#' @param chunksize The number of samples to take from the parameter space at each iteration
#' @param theta_sampler Function to take samples from the \eqn{Theta} parameter space. Default is \link{runif_dk_vects}. Must be a function of two parameters \code{d_k} a vector of dimensions, and \code{chunksize} the number of samples to take, and return a matrix with \code{sum(d_k)} columns and \code{chunksize} rows. See examples.
#' @param ga Logical, if TRUE, uses gradient ascent.
#' @param ga_gfactor Concentration parameter scale in the gradient ascent algorithm. A number or "adapt"
#' @param ga_lrate The gradient ascent learning rate
#' @param ga_restart_every Restart the gradient ascent after this number of iterations at a sample from \code{theta_sampler}
#' @param seed Seed for the random number generator. Can be set to NULL in which case no seed is set.
#'
#' @returns An object of class "htest", which is a list with the following elements:
#' \describe{
#' \item{estimate}{The value of the statistic at the observed data}
#' \item{p.value}{The p value}
#' \item{conf.int}{The upper and lower confidence limits}
#' \item{null.value}{The null hypothesis value provided by the user}
#' \item{alternative}{The type of test}
#' \item{method}{A description of the method}
#' \item{data.name}{The name of the data object provided by the user}
#' \item{p.sequence}{A list with two elements, p.null and p.alt containing the vector of p values at each iteration for the less than null and the greater than null. Used for assessing convergence. }
#' }
#'
#' @references Sachs, M.C., Gabriel, E.E. and Fay, M.P., 2024. Exact confidence intervals for functions of parameters in the k-sample multinomial problem. arXiv preprint arXiv:2406.19141.
#'
#'
#' @export
#' @examples
#' tau_ba <- function(theta) {
#'   theta1 <- theta[1:4]
#'   theta2 <- theta[5:8]
#'   sum(sqrt(theta1 * theta2))
#'   }
#' data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
#' xactonomial(data, tau_ba, psi_limits = c(0, 1), psi0 = .5,
#'   conf_int = FALSE, maxit = 15, chunksize = 200)
#'
#' # vectorized by row
#' tau_ba_v <- function(theta) {
#' theta1 <- theta[,1:4, drop = FALSE]
#' theta2 <- theta[,5:8, drop = FALSE]
#' rowSums(sqrt(theta1 * theta2))
#' }
#' data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
#' xactonomial(data, tau_ba_v, psi_limits = c(0, 1), psi0 = .5,
#'  conf_int = FALSE, maxit = 10, chunksize = 200)
#'
#'  # example of using theta_null_points
#'  # psi = 1/3 occurs when all probs = 1/3
#'  tau_max <- function(pp) {
#'    max(pp)
#'  }
#'
#' data <- list(c(13, 24, 13))
#'
#' xactonomial(data, tau_max, psi_limits = c(1 / 3, 1), psi0 = 1/ 3,
#'   conf_int = FALSE, theta_null_points = t(c(1/3, 1/3, 1/3)))
#'
#' ## in this case using p_value_limits improves confidence interval performance
#'
#'  xactonomial(data, tau_max, psi_limits = c(1 / 3, 1), psi0 = 1/ 3,
#'   conf_int = TRUE, theta_null_points = t(c(1/3, 1/3, 1/3)),
#'   p_value_limits = c(.1, 1e-8))
#'
#' ## specifying theta_sampler
#'
#' dirich_sampler <- function(d_k, chunksize){
#'  rdirich_dk_vects(chunksize, list(1:4 + 1))
#' }
#'
#' xactonomial(list(1:4),tau_max,
#'             psi_limits = c(0.25,1), psi0 = .5, conf_int = FALSE,
#'             theta_sampler = dirich_sampler)
#'

xactonomial <- function(data, f_param, statistic = NULL, psi0 = NULL,
                        alternative = c("two.sided", "less", "greater"),
                        psi_limits, theta_null_points = NULL, p_target = 1,
                        conf_int = TRUE, conf_level = .95, itp_maxit = 10, itp_eps = .005,
                        p_value_limits = NULL,
                        maxit = 50, chunksize = 500,
                        theta_sampler = runif_dk_vects,
                        ga = TRUE, ga_gfactor = "adapt", ga_lrate = .01,
                        ga_restart_every = 10,
                        seed = 503
                        ) {

  alternative <- match.arg(alternative)
  if(missing(psi_limits)) {
    stop("Please provide the lower and upper limits of the output of the f_param function as a vector in psi_limits")
  }

  if(!is.null(psi0)) {
    if(psi0 < min(psi_limits) | psi0 > max(psi_limits)) {
    stop("psi0 outside the possible range specified by psi_limits.")
  }}

  if(!is.null(psi0) & (any(abs(psi0 - psi_limits) < 1e-8) &
                       (is.null(theta_null_points) & isTRUE(all.equal(theta_sampler, runif_dk_vects))))) {


    warning("The psi0 value is on the boundary of the psi parameter space. The default Monte Carlo method may not perform well in this case. Consider specifying theta_null_points which is a matrix of theta vectors at which f_param(theta) = psi0, or specifying theta_sampler to sample directly from the null space.")

  }

  if(!is.null(seed)) {
    set.seed(seed)
  }

  alpha <- 1 - conf_level

  k <- length(data)
  d_k <- sapply(data, length)
  n_k <- sapply(data, sum)

  vtest <- runif_dk_vects(d_k, 5)
  f_is_vectorized <- tryCatch(length(f_param(vtest)) == 5, error = function(e) FALSE)

  sspace_size <- prod(sapply(data, \(dd) choose(length(dd) + sum(dd) - 1, length(dd) - 1)))

  SSpace <- lapply(data, \(x) matrix(sspace_multinom(length(x), sum(x)), ncol = length(x), byrow = TRUE))

  if(k == 1) {

    newX <- SSpace[[1]]
    sumX <- sum(newX[1,])
    logC <- lfactorial(sumX) - rowSums(apply(newX, 2, lfactorial))

    spacelist <- list(Sspace = newX, logC = logC)

  } else if(k == 2) {

    spacelist <- combinate(SSpace[[1]], SSpace[[2]])

  } else {

    spacelist <- combinate(SSpace[[1]], SSpace[[2]])
    for(i in 3:k) {

      spacelist <- combinate2(spacelist, SSpace[[i]])

    }

  }


  if(is.null(statistic)) { ## compute f_param with empirical proportions

    user_stat <- FALSE
    statistic <- function(df) {
      denom <- rep.int(n_k, d_k)
      if(is.matrix(df)) {
        pmat <- matrix(rep(denom, nrow(df)), nrow = nrow(df), ncol = length(denom), byrow = TRUE)
        if(f_is_vectorized) {
          f_param(df / pmat)
        } else {
          apply(df / pmat, 1, f_param)
        }
      } else {
        if(f_is_vectorized) {
          f_param(t(df / denom))
        } else {
          f_param(df / denom)
        }

      }

    }

  } else {
    user_stat <- TRUE
  }

  psi_obs <- statistic(unlist(data))

  psi_hat <- statistic(spacelist$Sspace)

  method <- "Monte Carlo multinomial test"
  pvalues <- if(!is.null(psi0)) {

    if(!is.null(theta_null_points) && any(abs(psi0 - psi_limits) < 1e-8)) {

      method <- "Exact multinomial test given a point null"
      II.lower <- psi_hat >= psi_obs
      II.upper <- psi_hat <= psi_obs
     c(max(calc_prob_null_fast(theta_null_points, spacelist$Sspace, spacelist$logC,  II.lower)),
       max(calc_prob_null_fast(theta_null_points, spacelist$Sspace, spacelist$logC,  II.upper)))


    } else {


    pvalue_psi0(psi0 = psi0, f_param = f_param, psi_hat = psi_hat, psi_obs = psi_obs,
                alternative = alternative,
                maxit = maxit, chunksize = chunksize,
                            p_target = p_target, SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                            d_k = d_k, f_is_vectorized = f_is_vectorized,
                            theta_sampler = theta_sampler,
                            ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                ga_restart_every = ga_restart_every)

    }

    } else NA


  confint <- if(conf_int) {
  flower <- function(x, f_param, psi_hat, psi_obs, maxit, chunksize, p_target, SSpacearr,
                     logC, d_k, f_is_vectorized,
                     theta_sampler, ga, ga_gfactor, ga_lrate, ga_restart_every){
    pest <- pvalue_psi0(x, f_param, psi_hat, psi_obs, maxit = maxit, chunksize = chunksize,
                p_target = p_target, SSpacearr = SSpacearr, logC = logC,
                d_k = d_k, f_is_vectorized = f_is_vectorized,
                theta_sampler = theta_sampler, alternative = "greater",
                ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                ga_restart_every = ga_restart_every, warn = FALSE)[1]
    (if(is.na(pest)) 0 else pest) - alpha / 2
    }

  if(!is.null(p_value_limits)){
    flower.boundary <- p_value_limits[1] - alpha / 2
  } else if(!is.null(theta_null_points)) {
    if(!is.null(psi0) & abs(psi0 - psi_limits[1]) < 1e-8) {

    II.lower <- psi_hat >= psi_obs

    flower.boundary <- max(calc_prob_null_fast(theta_null_points, spacelist$Sspace, spacelist$logC,  II.lower)) -
      alpha / 2
    } else {
      flower.boundary <- -alpha / 2
    }
  } else {
    flower.boundary <- -alpha / 2
  }


  fupper <- function(x, f_param, psi_hat, psi_obs, maxit, chunksize, p_target, SSpacearr, logC, d_k, f_is_vectorized,
                     theta_sampler, ga, ga_gfactor, ga_lrate, ga_restart_every){
    pest <- pvalue_psi0(x, f_param, psi_hat, psi_obs, maxit = maxit, chunksize = chunksize,
                p_target = p_target, SSpacearr = SSpacearr, logC = logC,
                d_k = d_k, f_is_vectorized = f_is_vectorized,
                theta_sampler = theta_sampler, alternative = "less",
                ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                ga_restart_every = ga_restart_every, warn = FALSE)[2]
    (if(is.na(pest)) 0 else pest) - alpha / 2
  }
  if(!is.null(p_value_limits)){
    fupper.boundary <- p_value_limits[2] - alpha / 2
  } else if(!is.null(theta_null_points)) {
    if(!is.null(psi0) & abs(psi0 - psi_limits[2]) < 1e-8) {
    II.upper <- psi_hat <= psi_obs
    fupper.boundary <- max(calc_prob_null_fast(theta_null_points, spacelist$Sspace, spacelist$logC,  II.upper)) -
      alpha / 2

    } else {
      fupper.boundary <- - alpha / 2
    }
  } else {
    fupper.boundary <- - alpha / 2
  }


  if(isTRUE(all.equal(psi_obs, min(psi_hat))) | flower.boundary > 0) {
    lower_limit <- psi_limits[1]
  } else {
    lower_limit <- itp_root(flower, psi_limits[1], psi_limits[2],
                            fa = flower.boundary, fb = 1 - alpha / 2, maxiter = itp_maxit,
                            eps = itp_eps,
                            f_param = f_param, psi_hat = psi_hat, psi_obs = psi_obs,
                            maxit = maxit, chunksize = chunksize,
                            p_target = alpha / 2 + 2 * itp_eps,
                            SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                            d_k = d_k, f_is_vectorized = f_is_vectorized,
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
                            eps = itp_eps,
                            f_param = f_param, psi_hat = psi_hat, psi_obs = psi_obs,
                            maxit = maxit, chunksize = chunksize,
                            p_target = alpha / 2 + 2 * itp_eps,
                            SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                            d_k = d_k, f_is_vectorized = f_is_vectorized,
                            theta_sampler = theta_sampler,
                            ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                            ga_restart_every = ga_restart_every)

    if((attr(upper_limit, "iter") %||%  itp_maxit) == itp_maxit){
      upper_limit <- psi_limits[2]
      warning("No root found for upper confidence limit, using upper boundary.")
    }

  }


  if(missing(p_value_limits) & is.null(theta_null_points)) {
  lowerp <-   pvalue_psi0(psi0 = lower_limit, f_param = f_param, psi_hat = psi_hat, psi_obs = psi_obs,
                         # alternative = "less",
                          maxit = maxit, chunksize = chunksize,
                          p_target = p_target, SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                          d_k = d_k, f_is_vectorized = f_is_vectorized,
                          theta_sampler = theta_sampler,
                          ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                          ga_restart_every = ga_restart_every, warn = FALSE)[1]


  if(lowerp > alpha/ 2) {
    warning("No root found for lower confidence limit, using lower boundary.")
    lower_limit <- psi_limits[1]

  }

  upperp <-   pvalue_psi0(psi0 = upper_limit, f_param = f_param, psi_hat = psi_hat, psi_obs = psi_obs,
                          #alternative = "greater",
                          maxit = maxit, chunksize = chunksize,
                          p_target = p_target, SSpacearr = spacelist$Sspace, logC = spacelist$logC,
                          d_k = d_k, f_is_vectorized = f_is_vectorized,
                          theta_sampler = theta_sampler,
                          ga = ga, ga_gfactor = ga_gfactor, ga_lrate = ga_lrate,
                          ga_restart_every = ga_restart_every, warn = FALSE)[2]


  if(upperp > alpha/ 2) {
    warning("No root found for upper confidence limit, using upper boundary.")
    upper_limit <- psi_limits[2]
  }
  }

  c(lower_limit, upper_limit)
  } else c(NA, NA)

  p.sequence <- attr(pvalues, "p.sequence")

  attr(confint, "conf.level") = 1 - alpha
  f_pname <- if(!user_stat) {
    deparse1(substitute(f_param))
  } else {
    deparse1(substitute(statistic))
  }
  names(psi_obs) <- f_pname

  res <- list(
    estimate = psi_obs,
    p.value = switch(alternative, "greater" = pvalues[1], "less" = pvalues[2],
                     "two.sided" = min(1, 2 * min(pvalues))),
    conf.int = confint,
    null.value = c(psi0 = psi0),
    alternative = alternative,
    method = method,
    data.name = deparse1(substitute(data)),
    p.sequence = p.sequence
  )
  class(res) <- "htest"
  res

}


#' Compute a p value for the test of psi <= psi0 and/or psi >= psi0
#'
#' @param psi0 The null hypothesis value for the parameter being tested.
#' @param f_param Function that takes in parameters and outputs a real
#'   valued number for each parameter. Can be vectorized rowwise for a matrix or not.
#' @param psi_hat The vector of psi values at each element of the sample space
#' @param psi_obs The observed estimate at the given data
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param maxit Maximum number of iterations of the Monte Carlo procedure
#' @param chunksize The number of samples to take from the parameter space at each iteration
#' @param p_target If a p-value is found that is greater than p_target, terminate the algorithm early.
#' @param SSpacearr The sample space matrix
#' @param logC The log multinomial coefficients for each row of the sample space
#' @param d_k The vector of dimensions
#' @param f_is_vectorized Is f_param vectorized by row?
#' @param theta_sampler Function to take samples from the \eqn{Theta} parameter space. Default is \link{runif_dk_vects}.
#' @param ga Logical, if TRUE, uses gradient ascent.
#' @param ga_gfactor Concentration parameter scale in the gradient ascent algorithm. A number or "adapt"
#' @param ga_lrate The gradient ascent learning rate
#' @param ga_restart_every Restart the gradient ascent after this number of iterations at a sample from
#' @param warn If TRUE, will give a warning if no samples from the null space are found
#' @returns A vector with two p-values, one for the lower, and one for the greater
#' @export
#' @examples
#'
#' sspace_3_5 <- matrix(sspace_multinom(3, 5), ncol = 3, byrow = TRUE)
#' f_max <- function(theta) max(theta)
#' logC <- apply(sspace_3_5, 1, log_multinom_coef, sumx = 5)
#' psi_hat <- apply(sspace_3_5, 1, \(x) f_max(x / sum(x)))
#' pvalue_psi0(.3, f_max, psi_hat, .4, maxit = 10, chunksize = 100,
#'  p_target = 1, SSpacearr = sspace_3_5, logC = logC, d_k = 3, warn = FALSE)

pvalue_psi0 <- function(psi0, f_param, psi_hat, psi_obs, alternative = "two.sided",
                        maxit, chunksize,
                        p_target,
                        SSpacearr, logC, d_k, f_is_vectorized = FALSE,
                        theta_sampler = runif_dk_vects,
                        ga = FALSE, ga_gfactor = 1, ga_lrate = .01,
                        ga_restart_every = 10, warn = TRUE
                        ) {

  already_warned <- !warn
  II.lower <- psi_hat >= psi_obs
  II.upper <- psi_hat <= psi_obs

  p.null <- rep(0, maxit)
  p.alt <- rep(0, maxit)

  theta_cands <- theta_sampler(d_k, chunksize)
  null_continue <- alt_continue <- TRUE
  if(alternative == "greater") alt_continue <- FALSE
  if(alternative == "less") null_continue <- FALSE
  null_stop <- alt_stop <- maxit
  never_null <- null_continue
  never_alt <- alt_continue
  for(i in 1:maxit) {

    if(isFALSE(null_continue) & isFALSE(alt_continue)) break

    this_theta <- theta_cands
    psi_theta <- if(f_is_vectorized) f_param(theta_cands) else apply(theta_cands, MARGIN = 1, f_param)
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
      probs_null <- calc_prob_null_fast(theta_null,
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
      probs_alt <- calc_prob_null_fast(theta_alt,
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
  if((isTRUE(never_null) | isTRUE(never_alt)) & !already_warned) {
    already_warned <- TRUE
    warning("Never found a value from the null space, results are not valid! Increase iterations, specify theta_null_points, or specify theta_sampler")
  }

  outnull <- if(isTRUE(never_null)) NA else max(p.null, na.rm = TRUE)
  outalt <- if(isTRUE(never_alt)) NA else max(p.alt, na.rm = TRUE)

  res <- c(null = outnull, alt = outalt)
  attr(res, "p.sequence") <- list(p.null = p.null[1:null_stop], p.alt = p.alt[1:alt_stop])
  res
}
