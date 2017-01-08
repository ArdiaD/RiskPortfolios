#' @name meanEstimation
#' @title Estimation of mean returns
#' @description Function which is used to compute the estimation of the mean returns.
#' @details The argument \code{control} is a list that can supply any of the following
#' components:
#' \itemize{ 
#' \item \code{type} method used to estimate the mean returns,
#' among \code{'naive'}, \code{'ewma'}, \code{'bs'} and \code{'mart'} where:
#' 
#' \code{'naive'} is used to compute the arithmetic mean of the returns.
#' 
#' \code{'ewma'} is used to compute the exponential weighted moving average
#' mean of the returns.  The data must be sorted from the oldest to the latest.
#' 
#' \code{'bs'} is used to compute the Bayes-Stein estimation.
#' 
#' \code{'mart'} is used to compute the Martinelli implied returns.
#' 
#' Default: \code{type = 'naive'}.
#' 
#' \item \code{lambda} decay parameter. Default: \code{lambda = 0.94}.
#' }
#' 
#' @param rets a \eqn{(T \times N)}{(T x N)} matrix of past returns.
#' @param control control parameters (see *Details*).
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of expected returns.
#' @author David Ardia, Kris Boudt and Jean-Philippe Gagnon Fleury.
#' @references 
#' Jorion, P. (2004).  Bayes-Stein Estimation for Portfolio
#' Analysis \emph{Journal of Finance and Quantitative Analysis} \bold{21}(3),
#' pp.279--292.
#' @keywords htest
#' @examples
#' # For the examples, we simply generate a 100 x 25 random matrix.
#' set.seed(3214)
#' T = 100
#' N = 25
#' rets = matrix(rnorm(T * N), nrow = T, ncol = N)
#' 
#' # Naive estimation of the mean
#' meanEstimation(rets)
#' 
#' # Naive estimation of the mean
#' meanEstimation(rets, control = list(type = 'naive'))
#' 
#' # Ewma estimation of the mean with default lambda = 0.94
#' meanEstimation(rets, control = list(type = 'ewma'))
#' 
#' # Ewma estimation of the mean with lambda = 0.9
#' meanEstimation(rets, control = list(type = 'ewma', lambda = 0.9))
#' 
#' # Martinelli's estimation of the mean
#' meanEstimation(rets, control = list(type = 'mart'))
#' 
#' # Bayes-Stein's estimation of the mean
#' meanEstimation(rets, control = list(type = 'bs'))
#' @export
meanEstimation <- function(rets, control = list()) {
  
  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x N) matrix")
  }
  
  ctr <- .ctrMean(control)
  
  if (ctr$type[1] == "naive") {
    mu <- .naiveMean(rets = rets)
  } else if (ctr$type[1] == "ewma") {
    mu <- .ewmaMean(rets = rets, lambda = ctr$lambda)
  } else if (ctr$type[1] == "bs") {
    mu <- .bsMean(rets = rets)
  } else if (ctr$type[1] == "mart") {
    mu <- .martMean(rets = rets)
  } else {
    stop("control$type is not well defined")
  }
  return(mu)
}

.ctrMean <- function(control = list()) {
  ## Function used to control the list input INPUTS control : a control
  ## ist The argument control is a list that can supply any of the
  ## following components type : default = 'naive' lambda : default = 0.94
  ## OUTPUTs control : a list
  if (!is.list(control)) {
    stop("control must be a list")
  }
  if (length(control) == 0) {
    control <- list(type = "naive", lambda = 0.94)
  }
  nam <- names(control)
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- c("naive", "ewma", "bs", "mart")
  }
  if (!("lambda" %in% nam) || is.null(control$lambda)) {
    control$lambda <- 0.94
  }
  return(control)
}

.naiveMean <- function(rets) {
  ## Compute the naive mean INPUTs rets : matrix (T x N) returns OUTPUTs
  ## mu : vector (N x 1) mean
  mu <- colMeans(rets)
  return(mu)
}

.ewmaMean <- function(rets, lambda) {
  ## Compute the exponential weighted moving average INPUTs rets : matrix
  ## (T x N) returns lambda : scalar OUTPUTs mu : vector (N x 1) mean NOTE
  ## the dates should be in ascending order, i.e. ordest at the top and
  ## newest at the bottom
  t <- dim(rets)[1]
  w <- lambda^(t:1)
  w <- w/sum(w)
  w <- matrix(data = w, nrow = t, ncol = dim(rets)[2], byrow = FALSE)
  mu <- colSums(w * rets)
  return(mu)
}

.bsMean <- function(rets) {
  ## Compute the Bayes-Stein mean INPUTs rets : matrix (T x N) returns
  ## OUTPUTs mu : vector (N x 1) mean !!! TOFIX !!! estime pour l'instant
  ## une matrice de var-cov de facon standard verifier avec david si on
  ## laisse le choix en input
  mu <- colMeans(rets)
  Sigma <- cov(rets)
  invSigma <- solve(Sigma)
  N <- dim(Sigma)[1]
  T <- length(mu)
  i <- rep(1, N)
  invSigmai <- crossprod(invSigma, i)
  w_min <- (invSigmai)/as.numeric(crossprod(i, invSigmai))
  mu_min <- crossprod(mu, w_min)
  invSigmaMu <- crossprod(invSigma, mu - mu_min)
  phi <- (N + 2)/((N + 2) + T * crossprod(mu - mu_min, invSigmaMu))
  phi <- max(min(phi, 1), 0)
  mu <- (1 - phi) * mu + phi * mu_min
  return(mu)
}

.martMean <- function(rets) {
  ## Compute the Martinelli mean INPUTs rets : matrix (T x N) returns
  ## OUTPUTs mu : vector (N x 1) mean
  mu <- apply(rets, 2, sd)
  return(mu)
}