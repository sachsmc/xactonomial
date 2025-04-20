test_that("confidence intervals for boundary problems", {


  psi_binom <- function(theta, d) {

    theta[1]

  }

  data_test <- list(c(0, 5))
  xact_bin <- xactonomial(data_test, psi_binom, psi_limits = c(0, 1), maxit = 100, chunksize = 50)

  expect_true(abs(xact_bin$conf.int[1] - 0) < 1e-8)


  psi_max <- function(pp) {

    max(pp)

  }

  data <- list(c(13, 24, 13))
  expect_warning(xactonomial(data, psi_max, psi_limits = c(1 / 3, 1), psi0 = 1/ 3, conf_int = FALSE, maxit = 10))
  run2 <- xactonomial(data, psi_max, psi_limits = c(1 / 3, 1), psi0 = 1/ 3,
              conf_int = FALSE, theta_null_points = t(c(1/3, 1/3, 1/3)))

  expect_true(run2$p.value > 1e-8)

  run3 <- xactonomial(data, psi_max, psi_limits = c(1 / 3, 1), psi0 = 1/ 2,
                      conf_int = TRUE, theta_null_points = t(c(1/3, 1/3, 1/3)))

  expect_false(abs(run3$conf.int[1] - 1/3) < 1e-8)

  run4 <- xactonomial(data, psi_max, psi_limits = c(1 / 3, 1),
                      conf_int = TRUE, p_value_limits = c(.1, 1e-8))

  expect_true(abs(run4$conf.int[1] - 1/3) < 1e-8)


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
  expect_lt(xactonomial(data, psi_bc, psi_limits = c(0, 1), psi0 = .5, conf_int = FALSE,
              maxit = 100, chunksize = 100)$p.value, .1)

  expect_gt(xactonomial(data, psi_bc, psi_limits = c(0, 1), psi0 = .5, alternative = "less", conf_int = FALSE,
              maxit = 100, chunksize = 50)$p.value, .1)

  expect_lt(xactonomial(data, psi_bc, psi_limits = c(0, 1), psi0 = .5, alternative = "greater", conf_int = FALSE,
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
  test <- xactonomial(data, psi_ba_v, psi_limits = c(0, 1), psi0 = .2, conf_int = FALSE,
                      maxit = 500, chunksize = 500,  ga = FALSE)
  expect_lt(test$p.value, .4)




})
