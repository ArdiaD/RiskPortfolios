context("Numerical values")

## The other test files check dimensions, summability and positivity only. The
## checks below compare the returned numbers against independently computed
## quantities: every one of them fails on version 2.1.7.

f.pRC <- function(w, Sigma) {
  Sw <- Sigma %*% w
  as.numeric(w * Sw)/as.numeric(crossprod(w, Sw))
}

test_that("erc equalises the risk contributions", {
  ## a covariance matrix with negative cross-correlations: an equity block and
  ## a bond block. The squared-deviation objective has a local minimum at
  ## "everything in the equity block", which is what a general-purpose
  ## optimizer started at the equally-weighted portfolio converges to.
  N = 10
  R = matrix(-0.4, N, N)
  R[1:5, 1:5] = 0.7
  R[6:10, 6:10] = 0.5
  diag(R) = 1
  s = c(rep(0.18, 5), rep(0.06, 5))
  Sigma = diag(s) %*% R %*% diag(s)

  for (cst in c('none', 'lo')) {
    w = optimalPortfolio(Sigma = Sigma, control = list(type = 'erc', constraint = cst))
    expect_equal(f.pRC(w, Sigma), rep(1/N, N), tolerance = 1e-8)
    expect_true(all(w > 0))
    expect_equal(sum(w), 1, tolerance = 1e-10)
  }

  ## and on a sample covariance matrix of the same shape
  set.seed(123)
  N = 20
  L = matrix(rnorm(N * N), N, N)/sqrt(N)
  Sigma = covEstimation(matrix(rnorm(250 * N), ncol = N) %*% L)
  w = optimalPortfolio(Sigma = Sigma, control = list(type = 'erc', constraint = 'lo'))
  expect_equal(f.pRC(w, Sigma), rep(1/N, N), tolerance = 1e-8)
})

test_that("the Bayes-Stein mean uses the number of observations", {
  ## 1 - phi -> 1 as T grows, so the estimator converges to the sample mean.
  set.seed(123)
  kept = sapply(c(100, 5000), function(T) {
    rets = matrix(rnorm(T * 10, mean = 0.005, sd = 0.05), ncol = 10)
    rets[, 1] = rets[, 1] + 0.02
    diff(range(meanEstimation(rets, control = list(type = 'bs'))))/
      diff(range(meanEstimation(rets)))
  })
  expect_true(kept[2] > 0.9)
  expect_true(kept[2] > kept[1])
})

test_that("maxdec uses the correlation matrix", {
  set.seed(123)
  rets = matrix(rnorm(250 * 10), ncol = 10) %*% diag(seq(0.5, 5, length.out = 10))
  Sigma = covEstimation(rets)
  w = optimalPortfolio(Sigma = Sigma, control = list(type = 'maxdec'))
  tmp = solve(cov2cor(Sigma), rep(1, 10))
  expect_equal(w, tmp/sum(tmp), tolerance = 1e-10)
  ## and is not the minimum-variance portfolio
  wmv = optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol'))
  expect_true(max(abs(w - wmv)) > 0.01)
})

test_that("gamma matters for the unconstrained mean-variance portfolio", {
  set.seed(123)
  rets = matrix(rnorm(250 * 10), ncol = 10)
  mu = meanEstimation(rets)
  Sigma = covEstimation(rets)
  w1 = optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv', gamma = 0.89))
  w2 = optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv', gamma = 20))
  expect_true(max(abs(w1 - w2)) > 1e-3)
  expect_equal(sum(w1), 1, tolerance = 1e-10)

  ## it must agree with the quadratic program without the bounds ...
  qp = quadprog::solve.QP(Dmat = 0.89 * Sigma, dvec = mu, Amat = matrix(1, 10, 1),
                          bvec = 1, meq = 1)$solution
  expect_equal(w1, qp, tolerance = 1e-8)
  ## ... and tend to the minimum-variance portfolio as gamma grows
  expect_equal(optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv', gamma = 1e8)),
               optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol')),
               tolerance = 1e-6)

  ## a portfolio with a negative 1'Sigma^-1 mu must not be sign-flipped
  w = optimalPortfolio(mu = -abs(mu), Sigma = Sigma, control = list(type = 'mv'))
  wmv = optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol'))
  expect_true(as.numeric(crossprod(w, -abs(mu))) > as.numeric(crossprod(wmv, -abs(mu))))
})

test_that("the ewma covariance weights sum to one", {
  set.seed(123)
  T = 20
  lambda = 0.94
  rets = matrix(rnorm(T * 3, sd = 0.02), ncol = 3)
  Sigma = covEstimation(rets, control = list(type = 'ewma', lambda = lambda))
  x = sweep(rets, 2, colMeans(rets), "-")
  ref = matrix(0, 3, 3)
  for (i in 1:T) {
    ref = (1 - lambda)/(1 - lambda^T) * outer(x[i, ], x[i, ]) + lambda * ref
  }
  expect_equal(Sigma, ref, tolerance = 1e-12)
})

test_that("'large' is 'lw'", {
  set.seed(123)
  rets = matrix(rnorm(250 * 10), ncol = 10)
  expect_equal(covEstimation(rets, control = list(type = 'large')),
               covEstimation(rets, control = list(type = 'lw')))

  ## a sample on which the untruncated intensity of the old 'large' turned
  ## negative, i.e. it shrank away from the prior: the two estimators differed
  ## by 1.6% there
  set.seed(211)
  r = matrix(rt(20 * 5, df = 3) * 0.02, ncol = 5)
  expect_equal(covEstimation(r, control = list(type = 'large')),
               covEstimation(r, control = list(type = 'lw')))
  expect_true(min(eigen(covEstimation(r, control = list(type = 'large')),
                        only.values = TRUE)$values) > 0)
})

test_that("riskeff works away from ten assets", {
  set.seed(123)
  for (N in c(3, 5, 9, 10, 25)) {
    rets = matrix(rnorm(200 * N), ncol = N)
    w = optimalPortfolio(Sigma = covEstimation(rets), semiDev = semidevEstimation(rets),
                         control = list(type = 'riskeff'))
    expect_equal(length(w), N)
    expect_equal(sum(w), 1, tolerance = 1e-6)
  }
  ## tied semideviations leave deciles empty
  set.seed(5)
  rets = matrix(rnorm(200 * 12), ncol = 12)
  rets[, 7:12] = rets[, 1:6]
  w = optimalPortfolio(Sigma = covEstimation(rets), semiDev = semidevEstimation(rets),
                       control = list(type = 'riskeff'))
  expect_false(anyNA(w))
})

test_that("the control arguments are validated", {
  set.seed(123)
  rets = matrix(rnorm(250 * 10), ncol = 10)
  mu = meanEstimation(rets)
  Sigma = covEstimation(rets)

  expect_error(covEstimation(rets, control = list(type = 'rtm')))
  expect_error(covEstimation(rets, control = list(type = 'ewma', lambda = 1)))
  expect_error(covEstimation(rets, control = list(type = 'ewma', lambda = -0.5)))
  expect_error(meanEstimation(rets, control = list(type = 'ewma', lambda = 0)))
  expect_error(covEstimation(rets, control = list(type = 'factor', K = 0)))
  expect_error(optimalPortfolio(mu = mu[1:3], Sigma = Sigma))
  expect_error(optimalPortfolio(Sigma = Sigma, semiDev = 1:3, control = list(type = 'riskeff')))
  expect_error(optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol',
                                constraint = 'gross', gross.c = 'abc')))
  expect_error(optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol',
                                constraint = 'gross', gross.c = 0.5)))
  expect_error(optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol',
                                constraint = 'user', LB = rep(0.2, 10), UB = rep(0.05, 10))))
  ## bounds incompatible with 1'w = 1, under any constraint
  expect_error(optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol',
                                constraint = 'user', LB = rep(0.2, 10), UB = rep(0.9, 10))))
  expect_error(optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol',
                                constraint = 'gross', LB = rep(0.2, 10), UB = rep(0.9, 10))))
  expect_error(optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol',
                                constraint = 'lo', UB = rep(0.05, 10))))
  ## bounds that bind on invvol must be honoured, not merely reported
  expect_warning(w <- optimalPortfolio(Sigma = Sigma, control = list(type = 'invvol',
                                       constraint = 'user', LB = rep(0.09, 10),
                                       UB = rep(0.11, 10))))
  expect_true(all(w >= 0.09 - 1e-9) && all(w <= 0.11 + 1e-9))
  expect_equal(sum(w), 1, tolerance = 1e-10)
  ## bounds that do not bind leave it untouched and silent
  expect_silent(w <- optimalPortfolio(Sigma = Sigma, control = list(type = 'invvol',
                                      constraint = 'user', LB = rep(0, 10), UB = rep(1, 10))))
  expect_equal(w, optimalPortfolio(Sigma = Sigma, control = list(type = 'invvol')))
})

test_that("the gross exposure constraint is respected", {
  set.seed(123)
  Sigma = covEstimation(matrix(rnorm(250 * 10), ncol = 10))
  for (gc in c(1.0, 1.2, 1.6, 3)) {
    w = optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol',
                         constraint = 'gross', gross.c = gc))
    expect_true(sum(abs(w)) <= gc + 1e-6)
    expect_equal(sum(w), 1, tolerance = 1e-6)
  }
})

test_that("the bounds are honoured on every path", {
  ## a feasible problem must not fail because the default starting value
  ## violates the bounds
  S = diag(4)
  w = optimalPortfolio(S, mu = 1:4, control = list(type = 'mv', constraint = 'gross',
                       LB = c(0.9, 0, 0, 0)))
  expect_true(all(w >= c(0.9, 0, 0, 0) - 1e-09))
  expect_equal(sum(w), 1, tolerance = 1e-08)
  w = optimalPortfolio(S, control = list(type = 'minvol', constraint = 'gross',
                       LB = c(0.9, 0, 0, 0)))
  expect_true(all(w >= c(0.9, 0, 0, 0) - 1e-09))

  ## clipping a solution into the box and then rescaling it can push a weight
  ## back out; the weights must satisfy both constraints on return
  ctr = RiskPortfolios:::.ctrPortfolio(2, list(type = 'minvol', constraint = 'user',
                                               LB = c(0.6, 0), UB = c(1, 1)))
  w = suppressWarnings(RiskPortfolios:::.finalizeWeights(c(0.6, 0.6), ctr))
  expect_true(all(w >= c(0.6, 0) - 1e-09))
  expect_equal(sum(w), 1, tolerance = 1e-10)

  ## .projectBox must be the Euclidean projection, not merely feasible
  set.seed(8)
  for (k in 1:50) {
    n = sample(2:12, 1)
    v = rnorm(n)
    lb = runif(n, -1, 0.5/n)
    ub = lb + runif(n, 0.05, 2/n)
    if (sum(lb) > 1 || sum(ub) < 1) next
    p = RiskPortfolios:::.projectBox(v, lb, ub)
    q = quadprog::solve.QP(diag(n), v, cbind(rep(1, n), diag(n), -diag(n)),
                           c(1, lb, -ub), meq = 1)$solution
    expect_equal(p, q, tolerance = 1e-07)
  }
})

test_that("a single asset does not break the estimators", {
  set.seed(123)
  rets = matrix(rnorm(100), ncol = 1)
  for (ty in c('naive', 'ewma', 'lw', 'const', 'cor', 'oneparm', 'diag', 'large', 'bs')) {
    Sigma = covEstimation(rets, control = list(type = ty))
    expect_equal(dim(Sigma), c(1L, 1L))
    expect_true(is.finite(Sigma[1, 1]) && Sigma[1, 1] > 0)
  }
  Sigma = covEstimation(rets)
  for (ty in c('minvol', 'invvol', 'erc', 'maxdiv', 'maxdec')) {
    expect_equal(optimalPortfolio(Sigma = Sigma, control = list(type = ty)), 1)
  }
})
