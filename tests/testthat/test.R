context("Cov dimensions") 

test_that("Cov dimensions", { 
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  Sigma = covEstimation(rets)
  expect_equal(dim(Sigma), c(25, 25)) 
}) 
