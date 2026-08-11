#' @name covEstimation
#' @aliases covEstimation
#' @title Covariance matrix estimation
#' @description Function which performs various estimations of covariance matrices.
#' @details The argument \code{control} is a list that can supply any of the following
#' components: 
#' \itemize{
#' \item \code{type} method used to compute the
#' covariance matrix, among \code{'naive'}, \code{'ewma'}, \code{'lw'},
#' \code{'factor'},\code{'const'}, \code{'cor'}, \code{'oneparm'},
#' \code{'diag'}, \code{'large'} and \code{'bs'} where:
#'
#' \code{'naive'} is used to compute
#' the naive (standard) covariance matrix.
#'
#' \code{'ewma'} is used to compute the exponential weighting moving average covariance matrix. The following formula is used
#' to compute the ewma covariance matrix:
#' \deqn{\Sigma := \frac{1-\lambda}{1-\lambda^T} \sum_{t=1}^{T} \lambda^{T-t} (r_t - \bar{r})(r_t - \bar{r})'}{Sigma
#' := (1-lambda)/(1-lambda^T) sum_t lambda^(T-t) (r[t] - rbar)(r[t] - rbar)'}
#' where \eqn{r_t} is the \eqn{(N \times 1)}{(N x 1)} vector of returns at time
#' \eqn{t}. The weights sum to one, so the estimator is the finite-sample
#' counterpart of the recursion \eqn{\Sigma_t := \lambda \Sigma_{t-1} + (1-\lambda) r_t r_t'}{Sigma[t]
#' := lambda * Sigma[t-1] + (1-lambda) r[t]r[t]'} started at zero.
#' Note that the data must be sorted from the oldest to the latest. See RiskMetrics (1996)
#'
#' \code{'factor'} is used to compute the covariance matrix estimation using a
#' K-factor approach. See Harman (1976).
#' 
#' \code{'lw'} is a weighted average of the sample covariance matrix and a
#' 'prior' or 'shrinkage target'. The prior is given by a one-factor model and
#' the factor is equal to the cross-sectional average of all the random
#' variables. See Ledoit and Wolf (2003).
#' 
#' \code{'const'} is a weighted average of the sample covariance matrix and a
#' 'prior' or 'shrinkage target'.  The prior is given by constant correlation
#' matrix. See Ledoit and Wolf (2002).
#' 
#' \code{'cor'} is a weighted average of the sample covariance matrix and a
#' 'prior' or 'shrinkage target'.  The prior is given by the constant
#' correlation covariance matrix given by Ledoit and Wolf (2003).
#' 
#' \code{'oneparm'} is a weighted average of the sample covariance matrix and a
#' 'prior' or 'shrinkage target'.  The prior is given by the one-parameter
#' matrix. All variances are the same and all covariances are zero. 
#' See Ledoit and Wolf (2004).
#' 
#' \code{'diag'} is a weighted average of the sample covariance matrix and a
#' 'prior' or 'shrinkage target'.  The prior is given by a diagonal matrix. 
#' See Ledoit and Wolf (2002).
#' 
#' \code{'large'} is an alias of \code{'lw'}, kept for backward compatibility.
#' The shrinkage intensity of the market-prior estimator is derived so as to
#' minimize the quadratic loss measured by the Frobenius norm, and is valid as
#' the number of variables and/or the number of observations go to infinity;
#' Monte-Carlo simulations show that it works well for values as low as 10. The
#' main advantage is that the estimator is guaranteed to be invertible and
#' well-conditioned even if variables outnumber observations.
#' See Ledoit and Wolf (2003, 2004).
#'
#' \code{'bs'} is the Bayes-Stein estimator for the covariance matrix given by
#' Jorion (1986).
#'
#' Default: \code{type = 'naive'}.
#'
#' \item \code{lambda} decay parameter, a single number in \eqn{(0, 1)}.
#' Default: \code{lambda = 0.94}.
#'
#' \item \code{K} number of factors to use when the K-factor approach is
#' chosen to estimate the covariance matrix. Default: \code{K = 1}.}
#'
#' @param rets a matrix \eqn{(T \times N)}{(T x N)} of returns.
#' @param control control parameters (see *Details*).
#' @return A \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @note Part of the code is adapted from the Matlab code by Ledoit and Wolf (2014).
#' @author David Ardia, Kris Boudt and Jean-Philippe Gagnon Fleury.
#' @references 
#' Jorion, P. (1986). 
#' Bayes-Stein estimation for portfolio analysis.
#' \emph{Journal of Financial and Quantitative Analysis} \bold{21}(3), pp.279-292. 
#' 
#' Harman, H.H. (1976)
#' \emph{Modern Factor Analysis}. 
#' 3rd Ed. Chicago: University of Chicago Press.
#' 
#' Ledoit, O., Wolf, M. (2002).  
#' Improved estimation of the covariance matrix of stock returns with an application to portfolio selection. 
#' \emph{Journal of Empirical Finance} \bold{10}(5), pp.603-621. 
#' 
#' Ledoit, O., Wolf, M. (2003).  
#' Honey, I Shrunk the Sample Covariance Matrix.
#' \emph{Journal of Portfolio Management} \bold{30}(4), pp.110-119. 
#' 
#' Ledoit, O., Wolf, M. (2004).  
#' A well-conditioned estimator for large-dimensional covariance matrices.
#' \emph{Journal of Multivariate Analysis} \bold{88}(2), pp.365-411. 
#' 
#' RiskMetrics (1996)
#' \emph{RiskMetrics Technical Document}.
#' J. P. Morgan/Reuters. 
#' @keywords htest
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets = Industry_10
#' 
#' # Naive covariance estimation
#' covEstimation(rets)
#' 
#' # Ewma estimation of the covariance with default lambda = 0.94
#' covEstimation(rets, control = list(type = 'ewma'))
#' 
#' # Ewma estimation of the covariance with default lambda = 0.90
#' covEstimation(rets, control = list(type = 'ewma', lambda = 0.9))
#' 
#' # Factor estimation of the covariance with dafault K = 1
#' covEstimation(rets, control = list(type = 'factor'))
#' 
#' # Factor estimation of the covariance with K = 3
#' covEstimation(rets, control = list(type = 'factor', K = 3))
#' 
#' # Ledot-Wolf's estimation of the covariance
#' covEstimation(rets, control = list(type = 'lw'))
#' 
#' # Shrinkage of the covariance matrix using constant correlation matrix
#' covEstimation(rets, control = list(type = 'const'))
#' 
#' # Shrinkage of the covariance matrix towards constant correlation matrix by
#' # Ledoit-Wolf.
#' covEstimation(rets, control = list(type = 'cor'))
#' 
#' # Shrinkage of the covariance matrix towards one-parameter matrix
#' covEstimation(rets, control = list(type = 'oneparm'))
#' 
#' # Shrinkage of the covariance matrix towards diagonal matrix
#' covEstimation(rets, control = list(type = 'diag'))
#'
#' # Shrinkage of the covariance matrix for large data set (alias of 'lw')
#' covEstimation(rets, control = list(type = 'large'))
#'
#' # Bayes-Stein estimation of the covariance
#' covEstimation(rets, control = list(type = 'bs'))
#' @export
#' @importFrom stats cor cov factanal median quantile sd
covEstimation <- function(rets, control = list()) {
  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x N) matrix")
  }
  
  ctr <- .ctrCov(control)
  
  if (ctr$type[1] == "naive") {
    Sigma <- .naiveCov(rets = rets)
  } else if (ctr$type[1] == "ewma") {
    Sigma <- .ewmaCov(rets = rets, lambda = ctr$lambda)
  } else if (ctr$type[1] == "lw") {
    Sigma <- .lwCov(rets = rets)
  } else if (ctr$type[1] == "factor") {
    Sigma <- .factorCov(rets = rets, K = ctr$K)
  } else if (ctr$type[1] == "const") {
    Sigma <- .constCov(rets = rets)
  } else if (ctr$type[1] == "cor") {
    Sigma <- .corCov(rets = rets)
  } else if (ctr$type[1] == "oneparm") {
    Sigma <- .oneparmCov(rets = rets)
  } else if (ctr$type[1] == "diag") {
    Sigma <- .diagCov(rets = rets)
  } else if (ctr$type[1] == "large") {
    Sigma <- .largeCov(rets = rets)
  } else if (ctr$type[1] == "bs") {
    Sigma <- .bsCov(rets = rets)
  } else {
    stop("control$type is not well defined")
  }
  return(Sigma)
}

.ctrCov <- function(control = list()) {
  #Function used to control the list input INPUTS control : a control
  #list The argument control is a list that can supply any of the
  #following components type : default = 'naive' lambda : default = 0.94
  #K : default = 1 OUTPUTs control : a list
  if (!is.list(control)) {
    stop("control must be a list")
  }
  nam <- names(control)
  type <- c("naive", "ewma", "lw", "factor", "const", "cor", "oneparm",
            "diag", "large", "bs")
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- type
  }
  if (!(control$type[1] %in% type)) {
    stop("control$type is not well defined")
  }
  if (!("lambda" %in% nam) || is.null(control$lambda)) {
    control$lambda <- 0.94
  }
  .checkLambda(control$lambda)
  if (!("K" %in% nam) || is.null(control$K)) {
    control$K <- 1
  }
  if (!is.numeric(control$K) || length(control$K) != 1L || control$K < 1 ||
      control$K != round(control$K)) {
    stop("'K' must be a single positive integer")
  }
  return(control)
}

.naiveCov <- function(rets) {
  ## Compute the naive covariance matrix INPUTs rets : matrix (T x N)
  ## returns OUTPUTs Sigma : matrix (N x N) covariance matrix DA here we
  ## could check that the dimension doesn't lead to bullshit outputs
  Sigma <- cov(rets)
  return(Sigma)
}

.ewmaCov <- function(rets, lambda) {
  ## Compute the exponential weighted moving average covariance matrix
  ## INPUTs rets : matrix (T x N) returns lambda : decay parameter
  ## OUTPUTs Sigma : matrix (N x N) covariance matrix
  ## Finite-sample ewma: the weights (1 - lambda)/(1 - lambda^t) * lambda^(t-i)
  ## sum to one, which is the recursion started at Sigma = 0. Starting it at
  ## cov(rets) instead would add a further lambda^t * cov(rets) on top of an
  ## already normalized estimator, and inflate the covariance by a factor
  ## (1 + lambda^t) -- 30% for t = 20 at the default lambda.
  t <- nrow(rets)
  mu <- colMeans(rets)
  shiftRets <- sweep(rets, 2, mu, "-")
  w <- (1 - lambda)/(1 - lambda^t) * lambda^((t - 1):0)
  Sigma <- crossprod(sqrt(w) * shiftRets)
  return(Sigma)
}

.constCov <- function(rets) {
  ## shrinkage of the covariance matrix using constant correlation matrix
  ## INPUTs rets : matrix (T x N) returns OUTPUTs Sigma : matrix (N x N)
  ## covariance matrix
  n <- dim(rets)[2]
  tmpMat <- matrix(rep(1, n^2), ncol = n)

  rho <- if (n > 1) mean(cor(rets)[lower.tri(tmpMat, diag = FALSE)]) else 0
  R <- rho * tmpMat
  diag(R) <- 1

  std <- apply(rets, 2, sd)
  diagStd <- diag(std, nrow = n)
  Sigma <- diagStd %*% R %*% diagStd

  return(Sigma)
}

.factorCov <- function(rets, K) {
  ## Compute the covariance matrix using K-factor approach INPUTs rets :
  ## matrix (T x N) returns K : [scalar] number of factors OUTPUTs Sigma :
  ## matrix (N x N) covariance matrix NOTE Matlab function is factoran
  ## (statistical toolbox)
  n <- dim(rets)[2]
  std <- apply(rets, 2, sd)

  fit <- factanal(rets, K)
  loading <- fit$loadings
  uniquenesses <- fit$uniquenesses

  R <- tcrossprod(loading) + diag(uniquenesses, nrow = n)
  diagStd <- diag(std, nrow = n)
  Sigma <- diagStd %*% R %*% diagStd

  return(Sigma)
}

.lwCovElement <- function(rets, type) {
  ## Computes the common elements of the functions for Ledoit-Wolf INPUTs
  ## rets : matrix (T x N) returns type : 'large', 'lw' or else OUTPUTs
  ## list : useful outputs
  t <- dim(rets)[1]
  n <- dim(rets)[2]
  mu <- colMeans(rets)
  shiftRets <- sweep(rets, 2, mu, "-")
  y <- shiftRets^2
  
  if (type == "large" || type == "lw") {
    # mkt = rowMeans(rets)
    mkt <- rowMeans(shiftRets)
    
    smple <- cov(cbind(rets, mkt)) * (t - 1)/t
    covmkt <- smple[1:n, n + 1]
    covmkt_ <- matrix(rep(covmkt, n), ncol = n, byrow = FALSE)
    varmkt <- as.numeric(smple[n + 1, n + 1])
    smple <- smple[-(n + 1), -(n + 1), drop = FALSE]
    
    prior <- outer(covmkt, covmkt)/varmkt
    diag(prior) <- diag(smple)
    
    lwCovElement <- list(t = t, n = n, mu = mu, shiftRets = shiftRets, 
                         y = y, smple = smple, mkt = mkt, covmkt = covmkt, covmkt_ = covmkt_, 
                         varmkt = varmkt, prior = prior)
  } else {
    smple <- (1/t) * crossprod(shiftRets)
    
    lwCovElement <- list(t = t, n = n, mu = mu, shiftRets = shiftRets, 
                         y = y, smple = smple, mkt = NULL, covmkt = NULL, covmkt_ = NULL, 
                         varmkt = NULL, prior = NULL)
  }
  return(lwCovElement)
}

.lwCov <- function(rets) {
  ## Shrinkage of the covariance matrix towards market INPUTs rets :
  ## matrix (T x N) returns OUTPUTs Sigma : matrix (N x N) covariance
  ## matrix
  ## Adapted from covMarket.m by Olivier Ledoit and Michael Wolf (2014)
  lwCovElement <- .lwCovElement(rets, type = "lw")
  
  t <- lwCovElement$t
  n <- lwCovElement$n
  mu <- lwCovElement$mu
  shiftRets <- lwCovElement$shiftRets
  mkt <- lwCovElement$mkt
  covmkt <- lwCovElement$covmkt
  varmkt <- lwCovElement$varmkt
  smple <- lwCovElement$smple
  prior <- lwCovElement$prior
  y <- lwCovElement$y
  z <- sweep(shiftRets, 1, mkt, "*")
  
  # Phi hat
  phiMat <- crossprod(y)/t - 2 * crossprod(shiftRets) * smple/t + smple^2
  phi <- sum(apply(phiMat, 2, sum))
  
  # Rho hat
  rhoMat1 <- 1/t * sweep(crossprod(y, z), 2, covmkt, "*")/varmkt
  rhoMat3 <- 1/t * crossprod(z) * outer(covmkt, covmkt)/varmkt^2
  rhoMat <- 2 * rhoMat1 - rhoMat3 - prior * smple
  diag(rhoMat) <- diag(phiMat)
  rho <- sum(apply(rhoMat, 2, sum))
  
  # Gamma hat
  gamma <- norm(smple - prior, "F")^2
  
  # shrinkage value
  shrinkage <- .shrinkage(phi - rho, gamma, t)
  
  # Sigma hat
  Sigma <- shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}

.largeCov <- function(rets) {
  ## INPUTs rets : matrix (T x N) returns OUTPUTs Sigma : matrix (N x N)
  ## covariance matrix
  ## 'large' used to carry its own transcription of covMarket.m in which every
  ## quantity was scaled by 1/(n t). The scaling cancels in the shrinkage
  ## intensity, so the estimator was numerically identical to 'lw' -- except
  ## that the intensity was not truncated to [0, 1] and could turn negative
  ## (anti-shrinkage) in small samples. It now delegates to .lwCov.
  return(.lwCov(rets))
}

.corCov <- function(rets) {
  ## Shrinkage of the covariance matrix towards constant correlation
  ## matrix by Ledoit-Wolf INPUTs rets : matrix (T x N) returns OUTPUTs
  ## Sigma : matrix (N x N) covariance matrix
  ## Adapted from covCor.m by Olivier Ledoit and Michael Wolf (2014)
  lwCovElement <- .lwCovElement(rets, type = "cor")
  
  t <- lwCovElement$t
  n <- lwCovElement$n
  mu <- lwCovElement$mu
  shiftRets <- lwCovElement$shiftRets
  y <- lwCovElement$y
  smple <- lwCovElement$smple
  var <- diag(smple)
  sqrtvar <- sqrt(var)
  outerSqrtVar <- outer(sqrtvar, sqrtvar)
  
  rBar <- if (n > 1) (sum(sum(smple/outerSqrtVar)) - n)/(n * (n - 1)) else 0
  prior <- rBar * outerSqrtVar
  diag(prior) <- var
  
  # phi hat
  phiMat <- crossprod(y)/t - 2 * crossprod(shiftRets) * smple/t + smple^2
  phi <- sum(apply(phiMat, 2, sum))
  
  # rho hat
  
  term1 <- crossprod(shiftRets^3, shiftRets)/t
  term2 <- sweep(smple, 1, diag(smple), "*")
  # term3 = sweep(smple, 1, var, '*') # lorsqu'on developpe on pourrait
  # annuler les term3 et term4 ... R term4 = sweep(smple, 1, var, '*')
  rhoMat <- term1 - term2  #- term3 + term4
  diag(rhoMat) <- 0
  rho <- sum(diag(phiMat)) + rBar * sum(sum(outer(1/sqrtvar, sqrtvar) * 
                                              rhoMat))
  
  # gamma hat
  gamma <- norm(smple - prior, type = "F")^2
  
  # shrinkage value
  shrinkage <- .shrinkage(phi - rho, gamma, t)
  
  # Sigma hat
  Sigma <- shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}


.diagCov <- function(rets) {
  ## Shrinks towards diagonal matrix INPUTs rets : matrix (T x N) returns
  ## OUTPUTs Sigma : matrix (N x N) covariance matrix
  lwCovElement <- .lwCovElement(rets, type = "diag")
  
  t <- lwCovElement$t
  n <- lwCovElement$n
  mu <- lwCovElement$mu
  shiftRets <- lwCovElement$shiftRets
  y <- lwCovElement$y
  smple <- lwCovElement$smple
  prior <- diag(diag(smple), nrow = n)
  
  # phi hat
  phiMat <- crossprod(y)/t - 2 * (crossprod(shiftRets)) * smple/t + smple^2
  phi <- sum(apply(phiMat, 2, sum))
  
  # rho hat
  rho <- sum(diag(phiMat))
  
  # gamma hat
  gamma <- norm(smple - prior, "F")^2
  
  # shrinkage value
  shrinkage <- .shrinkage(phi - rho, gamma, t)
  
  # Sigma hat
  Sigma <- shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}

.oneparmCov <- function(rets) {
  ## Shrinks towards one-parameter matrix INPUTs rets : matrix (T x N)
  ## returns OUTPUTs Sigma : matrix (N x N) covariance matrix
  ## Adapted from cov1para.m by Olivier Ledoit and Michael Wolf (2014)
  lwCovElement <- .lwCovElement(rets, type = "oneparm")
  
  t <- lwCovElement$t
  n <- lwCovElement$n
  mu <- lwCovElement$mu
  shiftRets <- lwCovElement$shiftRets
  smple <- lwCovElement$smple
  y <- lwCovElement$y
  meanvar <- mean(diag(smple))
  prior <- meanvar * diag(n)
  
  # phi hat
  phiMat <- crossprod(y)/t - 2 * (crossprod(shiftRets)) * smple/t + smple^2
  phi <- sum(apply(phiMat, 2, sum))
  
  # gamma hat
  gamma <- norm(smple - prior, type = "F")^2
  
  # shrinkage value
  shrinkage <- .shrinkage(phi, gamma, t)
  
  # Sigma hat
  Sigma <- shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}


.bsCov <- function(rets) {
  ## Compute the Bayes-Stein covariance matrix INPUTs rets : matrix (T x
  ## returns OUTPUTs Sigma : Matrix (N x N) covariance David : see
  ## Kolusheva, Daniela. (July 2008) Out-of-sample Performance of Asset
  ## Allocation Strategies
  t <- dim(rets)[1]
  n <- dim(rets)[2]

  mu <- colMeans(rets)
  Sigma <- cov(rets)
  invSigma <- solve(Sigma)

  i <- rep(1, n)
  invSigmai <- crossprod(invSigma, i)
  w_min <- (invSigmai)/as.numeric(crossprod(i, invSigmai))
  mu_min <- as.numeric(crossprod(mu, w_min))
  invSigmaMu <- crossprod(invSigma, mu - mu_min)
  phi <- (n + 2)/((n + 2) + t * as.numeric(crossprod(mu - mu_min, invSigmaMu)))
  phi <- max(min(phi, 1), 0)

  tau <- t * phi/(1 - phi)

  ## phi = 1 (all expected returns equal) gives tau = Inf; the limits of the
  ## two coefficients below are 1 and 1/t.
  cSigma <- if (is.finite(tau)) 1 + 1/(t + tau) else 1
  cPrior <- if (is.finite(tau)) tau/(t * (t + 1 + tau)) else 1/t

  Sigma <- Sigma * cSigma +
    cPrior * outer(i, i)/as.numeric(crossprod(i, invSigmai))

  return(Sigma)
}