#' @name impliedReturns
#' @title Implied returns estimation
#' @description Function which computes the implied returns.
#' @details The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{
#' \item \code{type} method used to compute the implied returns,
#' among \code{'regression'}, \code{'robust'}, \code{'constraint'} and \code{'rl'}
#' where:
#' 
#' \code{'regression'} is use to compute the implied expected returns. The
#' implied expected returns are solved by fitting the following regression:
#' 
#' \deqn{\tilde{\mu} = a + b[\hat{\Sigma}\hat{w}]}{\tilde{\mu} = a +
#' b[\hat{\Sigma}\hat{w}]}
#' 
#' The option \code{'reg'} can be used with \code{'regression'} if you want a
#' robust estimation or a constrained estimation. See \code{'reg'}.
#' 
#' \code{'bl'} is used to compute the Black-Litterman implied expected returns.
#' The Black-Litterman implied expected return is estimated by fitting the
#' following regression:
#' 
#' \deqn{\mu = l \iota + \gamma \Sigma w.}
#' 
#' \item{list('reg')}{is used when type = \code{'regression'}. \code{reg} can
#' supply either \code{'robust'} or \code{'constraint'}.
#' 
#' \code{'robust'} is used to compute the robust estimation of the implied
#' expected returns.
#' 
#' \code{'constraint'} is used to compute the implied expected returns with the
#' restriction that the coefficient is larger than zero.
#' 
#' The constraint implied expected return estimate \eqn{\tilde{\mu}} is given
#' by the following equation:
#' 
#' \deqn{\tilde{\mu} = a + b(\hat{\Sigma}\hat{w})}{\tilde{\mu} = a +
#' b(\hat{\Sigma}\hat{w})}} 
#' 
#' \item \code{gamma} is used for the computation of the Black-Litterman implied 
#' expected returns. Default: \code{gamma = 0.89}. 
#' }
#' 
#' @param rets a \eqn{(T \times N)}{(T x N)} matrix of returns.
#' @param mu a \eqn{(N \times 1)}{(N x 1)} vector of mean (expected returns).
#' @param Sigma a \eqn{(N \times N)}{(N x N)} matrix of covariances.
#' @param semiDev a \eqn{(N \times 1)}{(N x 1)} vector of semideviation.
#' @param w a \eqn{(N \times 1)}{(N x 1)} vector of portfolio weights.
#' @param control control parameters (see *Details*).
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of implied expected returns.
#' @author David Ardia <david.ardia@@unine.ch> and Jean-Philippe Gagnon Fleury.
#' @references Ardia, D., Boudt, K. (2015).  Implied expected returns and the
#' choice of a mean-variance efficient portfolio proxy.  \emph{Journal of
#' Portfolio Management} \bold{41} (4), pp.68--81.
#' 
#' Best, M. J., Grauer, R. R., (1985).  Capital asset pricing compatible with
#' observed market value weights.  \emph{Journal of Finance} \bold{40}(5),
#' pp.85--103.
#' 
#' Black, F., Litterman, R., (September-October 1992).  Global portfolio
#' optimization.  \emph{Financial Analyst Journal} \bold{48}(5), pp.28--43.
#' @keywords htest
#' @examples
#' # For the examples, we simply generate a 100 x 25 random matrix.
#' set.seed(3214)
#' T = 100
#' N = 25
#' rets = matrix(rnorm(T * N), nrow = T, ncol = N)
#' 
#' mu = meanEstimation(rets)
#' Sigma = covEstimation(rets)
#' w = rep(1, N)/N
#' 
#' # Computes the implied expected returns by Black-Litterman with gamma = 0.89.
#' impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = 'bl'))
#' 
#' # Computes the implied expected returns by Black-Litterman with gamma = 1.
#' impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = 'bl', gamma = 1))
#' 
#' # Computes the impled expected returns.
#' impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = 'regression'))
#' 
#' # Compute the robust implied expected returns. impliedReturns(mu = mu, Sigma =
#' # Sigma, w = w, control = list(type = 'regression', reg = 'robust'))
#' 
#' # Compute the constraint implied expected returns. impliedReturns(mu = mu, Sigma
#' # = Sigma, w = w, control = list(type = 'regression', reg = 'constraint'))
#' @export
impliedReturns <- function(rets = NULL, mu = NULL, Sigma = NULL, semiDev = NULL, 
                           w = NULL, control = list()) {
  
  ctr <- .ctrImpliedReturns(control)
  
  if (is.null(mu)) {
    mu <- meanEstimation(rets = rets, control = ctr$control.mean)
  }
  if (is.null(Sigma) & is.null(semiDev)) {
    Sigma <- covEstimation(rets = rets, control = ctr$control.cov)
  }
  if (is.null(semiDev) & is.null(Sigma)) {
    semiDev <- semidevEstimation(rets, control = ctr$control.semidev)
  }
  if (is.null(w)) {
    w <- optimalPortfolio(mu = mu, Sigma = Sigma, semiDev = semiDev, 
                          control = ctr$control.ptf)
  }
  
  if (ctr$type[1] == "bl") {
    imu <- .blImpliedReturns(Sigma = Sigma, w = w, gamma = ctr$gamma[1])
  } else if (ctr$type[1] == "regression") {
    if (ctr$reg[1] == "standard") {
      imu <- .regImpliedReturns(mu = mu, Sigma = Sigma, w = w)
    } else if (ctr$reg[1] == "robust") {
      imu <- .robregImpliedReturns(mu = mu, Sigma = Sigma, w = w)
    } else if (ctr$reg[1] == "constraint") {
      imu <- .constraintImpliedReturns(mu = mu, Sigma = Sigma, w = w)
    } else {
      stop("control$reg not well defined")
    }
  } else {
    stop("control$type not well defined")
  }
  
  out <- list(imu = imu, mu = mu, Sigma = Sigma, w = w)
  return(out)
}

.ctrImpliedReturns <- function(control = list()) {
  
  if (!is.list(control)) {
    stop("control must be a list")
  }
  if (length(control) == 0) {
    control <- list(type = "regression", reg = "standard", control.mean = NULL, 
                    control.cov = NULL, control.ptf = NULL)
  }
  nam <- names(control)
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- c("regression", "bl")
  }
  if (!("gamma" %in% nam) || is.null(control$gamma)) {
    control$gamma <- .ctrPortfolio()$gamma
  }
  if (!("reg" %in% nam) || is.null(control$reg)) {
    control$reg <- c("standard", "robust", "constraint")
  }
  if (!("control.mean" %in% nam) || is.null(control$control.mean)) {
    control$control.mean <- .ctrMean()
  }
  if (!("control.cov" %in% nam) || is.null(control$control.cov)) {
    control$control.cov <- .ctrCov()
  }
  if (!("control.semidev" %in% nam) || is.null(control$control.semidev)) {
    control$control.semidev <- .ctrSemidev()
  }
  if (!("control.ptf" %in% nam) || is.null(control$control.ptf)) {
    control$control.ptf <- .ctrPortfolio()
  }
  return(control)
}

.blImpliedReturns <- function(Sigma, w, gamma) {
  
  imu <- gamma * crossprod(Sigma, w)
  return(imu)
}

.regImpliedReturns <- function(mu, Sigma, w) {
  
  # DA this is unconstrained for the moment
  imu <- as.vector(lm(mu ~ (Sigma * w))$fitted)
  return(imu)
}

.robregImpliedReturns <- function(mu, Sigma, w) {
  
  # DA robust regression here
  imu <- as.vector(MASS::rlm(mu ~ (Sigma * w))$fitted)
  return(imu)
}

.constraintImpliedReturns <- function(mu, Sigma, w) {
  
  x <- crossprod(Sigma, w)
  xBar <- mean(x)
  shiftX <- x - xBar
  yBar <- mean(mu)
  shiftY <- mu - yBar
  b <- max(crossprod(shiftX, shiftY)/sum(x^2), 0)
  a <- yBar - b * xBar
  muT <- a + b * x
  return(muT)
}
