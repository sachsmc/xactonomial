test_that("theta sampling functions work", {

  psi_max <- function(pp) { do.call(pmax.int, lapply(1:ncol(pp), \(i) pp[,i])) }

  psi_is_vectorized <- TRUE
  true_psi_max <- psi_max(t(c(.4, .4, .2)))

  sample_data2 <- function(n) {

    list(rmultinom(1, n, prob = c(.4, .4, .2))[, 1])

  }

  data <- list(c(11, 7, 2))

  k <- length(data)
  d_k <- sapply(data, length)

  tmpdat <- lapply(data, \(x) x / sum(x)) |> unlist()
  psi_obs <- psi_max(matrix(tmpdat, nrow = 1))

  SSpace <- lapply(data, \(x) sspace_multinom(length(x), sum(x)))

  newX <- SSpace[[1]]
  sumX <- sum(newX[1,])
  logC <- lfactorial(sumX) - rowSums(apply(newX, 2, lfactorial))

  spacelist <- list(Sspace = newX, Sprobs = newX / sumX, logC = logC)

  psi_hat <- psi_max(spacelist$Sprobs)

  set.seed(401)
  niter <- 0
  pees <- NULL
  p.unif <- 0
  repeat{
    p.unif <- max(p.unif, pvalue_psi0(psi0 = .4, psi = psi_max, psi_hat = psi_hat,
                          psi_obs = psi_obs, maxit = 1, chunksize = 200,
                          lower = TRUE, target = 1,
                          SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                          psi_v = TRUE))
    niter <- niter + 1
    pees <- c(pees, p.unif)
    if(0.255 - p.unif < .001) break

  }


  set.seed(401)
  niterd <- 0
  peesd <- NULL
  p.unifd <- 0
  dirich_param <- list(data[[1]] + c(sqrt(20)))
  repeat{
    p.unifd <- max(p.unifd, pvalue_psi0(psi0 = .4, psi = psi_max, psi_hat = psi_hat,
                                      psi_obs = psi_obs, maxit = 1, chunksize = 200,
                                      lower = TRUE, target = 1,
                                      SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                                      psi_v = TRUE, sample_theta = \(d_k, cs) rdirich_dk_vects(d_k, cs, dirich_param)))
    niterd <- niterd + 1
    peesd <- c(peesd, p.unifd)
    if(0.255 - p.unifd < .001) break

  }

  expect_true(abs(p.unifd - p.unif) < .001)
  expect_true(niter > niterd)

})
