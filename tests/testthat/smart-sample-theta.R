test_that("theta sampling functions work", {

  psi_max <- function(pp) { do.call(pmax.int, lapply(1:ncol(pp), \(i) pp[,i])) }

  psi_is_vectorized <- TRUE
  true_psi_max <- psi_max(t(c(.4, .4, .2)))

  sample_data2 <- function(n) {

    list(rmultinom(1, n, prob = c(.4, .4, .2))[, 1])

  }


  data <- list(c(11, 7, 2))
  #data <- sample_data(c(10,5))

  k <- length(data)
  d_k <- sapply(data, length)

  tmpdat <- lapply(data, \(x) x / sum(x)) |> unlist()
  psi_obs <- psi_max(matrix(tmpdat, nrow = 1))

  SSpace <- lapply(data, \(x) sspace_multinom(length(x), sum(x)))
  #spacelist <- combinate(SSpace[[1]], SSpace[[2]])
  #SSpace <- spacelist$Sspace

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
                          psi_obs = psi_obs, maxit = 1, chunksize = 2000,
                          lower = TRUE, target = 1,
                          SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                          psi_v = TRUE))
    niter <- niter + 1
    pees <- c(pees, p.unif)
    if(0.2555 - p.unif < .0005) break
    if(niter > 1e6) break
  }


  set.seed(401)
  niterd <- 0
  peesd <- NULL
  p.unifd <- 0
  dirich_param <- list(ceiling(data[[1]] + c(sqrt(20))))
  repeat{
    p.unifd <- max(p.unifd, pvalue_psi0(psi0 = .4, psi = psi_max, psi_hat = psi_hat,
                                      psi_obs = psi_obs, maxit = 1, chunksize = 2000,
                                      lower = TRUE, target = 1,
                                      SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                                      psi_v = TRUE, sample_theta = \(d_k, cs) rdirich_dk_vects(cs, dirich_param)))
    niterd <- niterd + 1
    peesd <- c(peesd, p.unifd)
    if(0.2555 - p.unifd < .0005) break

  }

  set.seed(401)
  niterg <- 0
  peesg <- NULL
  p.unifg <- 0

  repeat{
    p.unifg <- max(p.unifg, pvalue_psi0(psi0 = .4, psi = psi_max, psi_hat = psi_hat,
                                        psi_obs = psi_obs, maxit = 1, chunksize = 2000,
                                        lower = TRUE, target = 1,
                                        SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                                        psi_v = TRUE, gradient_ascent = TRUE, gamma = 10))
    niterg <- niterg + 1
    peesg <- c(peesg, p.unifg)
    if(0.2555 - p.unifg < .0005) break

  }

  expect_true(abs(p.unifd - p.unif) < .001)
  expect_true(abs(p.unifg - p.unif) < .001)
  expect_true(niter > niterd)
  expect_true(niter > niterg)

  # plot(pees, type = "l")
  # lines(peesd)
  # lines(peesg)

})


test_that("reviewer 1 problem is solved", {

  psi <- function(theta) {

    theta[, 1]

  }

  data <- list(c(20, rep(0, 9)))

  k <- length(data)
  d_k <- sapply(data, length)

  tmpdat <- lapply(data, \(x) x / sum(x)) |> unlist()
  psi_obs <- psi(matrix(tmpdat, nrow = 1))

  SSpace <- lapply(data, \(x) sspace_multinom(length(x), sum(x)))

  newX <- SSpace[[1]]
  sumX <- sum(newX[1,])
  logC <- lfactorial(sumX) - rowSums(apply(newX, 2, lfactorial))

  spacelist <- list(Sspace = newX, Sprobs = newX / sumX, logC = logC)

  psi_hat <- psi(spacelist$Sprobs)

  set.seed(401)
  niter <- 0
  pees <- NULL
  p.unif <- 0
  start <- Sys.time()
  repeat{
    p.unif <- max(p.unif, pvalue_psi0(psi0 = .1 ^(1 / 20), psi = psi, psi_hat = psi_hat,
                                      psi_obs = psi_obs, maxit = 1, chunksize = 2000,
                                      lower = TRUE, target = 1,
                                      SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                                      psi_v = TRUE))
    niter <- niter + 1
    pees <- c(pees, p.unif)
    if(0.1 - p.unif < .001) break
    if(niter > 1e6) break
  }
  secs <- Sys.time() - start

  set.seed(401)
  niterd <- 0
  peesd <- 1e-4
  p.unifd <- 0
  dirich_param <- list(ceiling(data[[1]] + .5 * (sum(data[[1]]))))
  start <- Sys.time()
  repeat{
    p.unifd <- max(p.unifd, pvalue_psi0(psi0 = .1^(1 / 20), psi = psi, psi_hat = psi_hat,
                                        psi_obs = psi_obs, maxit = 1, chunksize = 500, target = 1,
                                        SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                                        psi_v = TRUE, sample_theta = \(d_k, cs) rdirich_dk_vects(cs, dirich_param))[1])
    niterd <- niterd + 1
    peesd <- c(peesd, p.unifd)
    if(0.1 - p.unifd < .001) break

  }
  secsd <- Sys.time() - start

  set.seed(401)
  niterg <- 0
  peesg <- NULL
  p.unifg <- 1e-4
  dirich_param <- list(ceiling(data[[1]] + 1))
  start <- Sys.time()
  repeat{
    p.unifg <- max(p.unifg, pvalue_psi0(psi0 = .1^(1 / 20), psi = psi, psi_hat = psi_hat,
                                        psi_obs = psi_obs, maxit = 5, chunksize = 2000,
                                        target = 1,
                                        SSpacearr = spacelist$Sspace, logC = spacelist$logC, d_k = d_k,
                                        psi_v = TRUE, gradient_ascent = TRUE, gamma = 10,
                                        sample_theta = \(d_k, cs) rdirich_dk_vects(cs, dirich_param)))
    niterg <- niterg + 100
    peesg <- c(peesg, p.unifg)
    if(0.1 - p.unifg < .001) break

  }
  secsg <- Sys.time() - start

})


