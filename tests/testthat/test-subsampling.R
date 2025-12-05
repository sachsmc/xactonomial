test_that("Subsampling the sample space", {


  data1 <- list(A = c(40, 54, 61))
  fullssp <- matrix(sspace_multinom(3, sum(data1$A)), ncol = 3, byrow = TRUE)
  sampssp <- matrix(sspace_multinom_sample(3, sum(data1$A), 50), ncol = 3, byrow = TRUE)

  system.time(sspace_multinom_sample(3, sum(data1$A), 50))
  system.time(sspace_multinom(3, sum(data1$A)))

  expect_true(nrow(sampssp) == 51)

  datasm <- list(Zarr = c(21, 0, 1, 10),
               Zadv = c(15, 19, 3, 3))

  sspace_size <- prod(sapply(datasm, \(dd) choose(length(dd) + sum(dd) - 1, length(dd) - 1)))
  SSpace <- lapply(datasm, \(x) matrix(sspace_multinom(length(x), sum(x)), ncol = length(x), byrow = TRUE))


  databg <- list(Zarr = c(81, 0, 1, 10, 0, 0),
                 Zadv = c(15, 69, 3, 3, 15, 3),
                 Zsep = c(21, 4, 62, 5, 1, 20))



})
