test_that("confidence intervals for boundary problems", {


  psi_binom <- function(theta, d) {

    theta[1]

  }

  data_test <- list(c(0, 5))
  xact_bin <- xactonomial(psi_binom, data_test, psi_limits = c(0, 1), maxit = 100, chunksize = 50)

  expect_equal(xact_bin$conf.int[1], 0)

})


test_that("two sided tests", {



  psi_bc <- function(theta, d = 2) {

    gpsz <- length(theta) / d
    cprod <- rep(1, gpsz)
    j <- 1
    for(i in 1:d) {
      cprod <- cprod * theta[j:(j + gpsz - 1)]
      j <- j + gpsz
    }

    sum((cprod)^(1/d))

  }

  data <- list(T1 = c(2,1,2,1), T2 = c(0,1,3,3))
  set.seed(1)
  expect_lt(xactonomial(psi_bc, data, psi_limits = c(0, 1), psi0 = .5,
              maxit = 400, chunksize = 50)$p.value, .1)

  expect_gt(xactonomial(psi_bc, data, psi_limits = c(0, 1), psi0 = .5, alternative = "less",
              maxit = 100, chunksize = 20)$p.value, .1)

  expect_lt(xactonomial(psi_bc, data, psi_limits = c(0, 1), psi0 = .5, alternative = "greater",
              maxit = 400, chunksize = 50)$p.value, .1)

})


test_that("three samples", {
  psi_ba_v <- function(theta) {
  theta1 <- theta[,1:3, drop = FALSE]
  theta2 <- theta[,4:6, drop = FALSE]
  theta3 <- theta[,7:9, drop = FALSE]
  rowSums((theta1 * theta2 * theta3)^(1/3))
  }
  data <- list(T1 = c(2,1,1), T2 = c(0,1,3), T3 = c(1, 3, 0))

  set.seed(5)
  expect_lt(xactonomial(psi_ba_v, data, psi_limits = c(0, 1), psi0 = .1, conf.int = FALSE,
              maxit = 50, chunksize = 20, psi_is_vectorized = TRUE)$p.value, .4)




})
