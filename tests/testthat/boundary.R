test_that("confidence intervals for boundary problems", {


  psi_binom <- function(theta, d) {

    theta[1]

  }

  data_test <- list(c(0, 5))
  xact_bin <- xactonomial(psi_binom, data_test, psi_limits = c(0, 1), maxit = 100, chunksize = 50)

  expect_equal(xact_bin$conf.int[1], 0)

})
