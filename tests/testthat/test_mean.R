context("Mean") 

test_that("Mean dimensions", { 
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  
  mu = meanEstimation(rets)
  expect_equal(length(mu), N)
  
  mu = meanEstimation(rets, control = list(type = 'naive'))
  expect_equal(length(mu), N)
  
  mu = meanEstimation(rets, control = list(type = 'ewma'))
  expect_equal(length(mu), N)
  
  mu = meanEstimation(rets, control = list(type = 'ewma', lambda = 0.9))
  expect_equal(length(mu), N)
  
  mu = meanEstimation(rets, control = list(type = 'mart'))
  expect_equal(length(mu), N)
  
  mu = meanEstimation(rets, control = list(type = 'bs'))
  expect_equal(length(mu), N)
})