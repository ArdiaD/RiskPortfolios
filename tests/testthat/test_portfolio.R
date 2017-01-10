context("Portfolios") 

test_that("Portfolios", { 
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  
  f.test.pos = function(W) {
    out = all(w >= 0)
    return(out)
  }
  
  f.test.sum = function(w) {
    out = abs(sum(w) - 1) < 1e-4
    return(out)
  }
  
  mu = meanEstimation(rets)
  Sigma = covEstimation(rets)
  semiDev = semidevEstimation(rets)
  LB = rep(0.02, N)
  UB = rep(0.95, N)
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma)
  expect_equal(length(w), N)
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(gamma = 1))
  expect_equal(length(w), N)
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(type = 'mv'))
  expect_equal(length(w), N)
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(type = 'mv', constraint = 'none'))
  expect_equal(length(w), N)
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(type = 'mv', constraint = 'lo'))
  expect_equal(length(w), N)
  expect_true(f.test.pos(w))
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(type = 'mv', constraint = 'user', LB = LB, UB = UB))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(type = 'mv', constraint = 'gross'))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(type = 'mv', constraint = 'gross', LB = LB, UB = UB))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(mu = mu, Sigma = Sigma, 
                       control = list(type = 'mv', constraint = 'gross', gross.c = 1.2))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'minvol'))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'minvol', constraint = 'none'))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'minvol', constraint = 'lo'))
  expect_equal(length(w), N)
  expect_true(f.test.pos(w))
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'minvol', constraint = 'user', LB = LB, UB = UB))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'minvol', constraint = 'gross'))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'minvol', constraint = 'gross', gross.c = 1.2))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'invvol'))
  expect_equal(length(w), N)
  expect_true(f.test.pos(w))
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'erc'))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'erc', constraint = 'lo'))
  expect_equal(length(w), N)
  expect_true(f.test.pos(w))
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'erc', constraint = 'user', LB = LB, UB = UB))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'maxdiv'))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'maxdiv', constraint = 'lo'))
  expect_equal(length(w), N)
  expect_true(f.test.pos(w))
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, 
                       control = list(type = 'maxdiv', constraint = 'user', LB = LB, UB = UB))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, semiDev = semiDev, 
                       control = list(type = 'riskeff'))
  expect_equal(length(w), N)
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, semiDev = semiDev, 
                       control = list(type = 'riskeff', constraint = 'lo'))
  expect_equal(length(w), N)
  expect_true(f.test.pos(w))
  expect_true(f.test.sum(w))
  
  w = optimalPortfolio(Sigma = Sigma, semiDev = semiDev, 
                       control = list(type = 'riskeff', constraint = 'user', LB = LB, UB = UB))
  expect_equal(length(w), N)
  expect_true(f.test.pos(w))
  expect_true(f.test.sum(w))
})