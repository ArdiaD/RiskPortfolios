context("Covariance") 

test_that("Cov dimensions", { 
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  
  Sigma = covEstimation(rets)
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "ewma"))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "ewma", lambda = 0.9))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "factor"))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "factor", K = 3))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "lw"))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "const"))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "cor"))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "oneparm"))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "diag"))
  expect_equal(dim(Sigma), c(N, N))
  
  Sigma = covEstimation(rets, control = list(type = "large"))
  expect_equal(dim(Sigma), c(N, N))
}) 

test_that("Cov symmetry", { 
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  
  Sigma = covEstimation(rets)
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "ewma"))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "ewma", lambda = 0.9))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "factor"))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "factor", K = 3))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "lw"))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "const"))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "cor"))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "oneparm"))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "diag"))
  expect_true(isSymmetric(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "large"))
  expect_true(isSymmetric(Sigma))
}) 

test_that("Cov positive definite", { 
  
  f.test.pos = function(X) {
    out = min(eigen(X)$values) > 0
    return(out)
  }
  
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  
  Sigma = covEstimation(rets)
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "ewma"))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "ewma", lambda = 0.9))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "factor"))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "factor", K = 3))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "lw"))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "const"))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "cor"))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "oneparm"))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "diag"))
  expect_true(f.test.pos(Sigma))
  
  Sigma = covEstimation(rets, control = list(type = "large"))
  expect_true(f.test.pos(Sigma))
}) 
