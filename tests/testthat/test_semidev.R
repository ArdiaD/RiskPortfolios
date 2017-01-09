context("Semideviation") 

test_that("Dimensions", { 
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  
  sig = semidevEstimation(rets)
  expect_equal(length(sig), N)
  
  sig = semidevEstimation(rets, control = list(type = 'naive'))
  expect_equal(length(sig), N)
  
  sig = semidevEstimation(rets, control = list(type = 'ewma'))
  expect_equal(length(sig), N)
  
  sig = semidevEstimation(rets, control = list(type = 'ewma', lambda = 0.9))
  expect_equal(length(sig), N)
}) 

test_that("Positivity", { 
  
  f.test.pos = function(x) {
    out = all(x >= 0)
    return(out)
  }
  
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  
  sig = semidevEstimation(rets)
  expect_true(f.test.pos(sig))
  
  sig = semidevEstimation(rets, control = list(type = 'naive'))
  expect_true(f.test.pos(sig))
  
  sig = semidevEstimation(rets, control = list(type = 'ewma'))
  expect_true(f.test.pos(sig))
  
  sig = semidevEstimation(rets, control = list(type = 'ewma', lambda = 0.9))
  expect_true(f.test.pos(sig))
  
  Sigma = covEstimation(rets)
  expect_true(f.test.pos(sig))
}) 
