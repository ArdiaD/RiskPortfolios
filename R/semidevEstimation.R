#' @name semidevEstimation
#' @title Estimation of the semideviation
#' @description Function which computes the semideviation.
#' @details The argument \code{control} is a list that can supply any of the following
#' components:
#' 
#' \itemize{
#' \item \code{type} method used to compute the semideviation
#' vector, among \code{'naive'} and \code{'ewma'} where:
#' 
#' \code{'naive'} is used to compute the simple semideviation.
#' 
#' \code{'ewma'} is used to compute the exponential weighted moving average
#' semideviation. The data must be sorted from the oldest to the latest.
#' 
#' The semideviation for one stock is computed as follows. First we select the
#' returns which are smaller than the average of the past returns; we get a new
#' vector of dimension \eqn{K \times 1, K \le N}. Then, the weight \eqn{w_i}
#' for each observation at its corresponding time \eqn{t} is computed as \eqn{w
#' = \lambda^{t}}. We obtain a \eqn{K \times 1}{Kx1} vector. The vector of
#' weights is then normalized.  Finally, the semideviation is obtained as the
#' weighted standard deviation.
#' 
#' Default: \code{type = 'naive'}.
#' 
#' \item \code{lambda} decay parameter. Default: \code{lambda = 0.94}.
#' }
#' 
#' @param rets a \eqn{(T \times N)}{(T x N)} matrix of past returns.
#' @param control control parameters (see *Details*).
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of semideviations.
#' @author David Ardia, Kris Boudt and Jean-Philippe Gagnon Fleury.
#' @keywords htest
#' @examples
#' set.seed(3214)
#' T = 100
#' N = 25
#' rets = matrix(rnorm(T * N), nrow = T, ncol = N)
#' 
#' # Naive semideviation estimation
#' semidevEstimation(rets)
#' 
#' # Naive estimation of the semideviation
#' semidevEstimation(rets, control = list(type = 'naive'))
#' 
#' # Ewma estimation of the semideviation. Default lambda = 0.94
#' semidevEstimation(rets, control = list(type = 'ewma'))
#' 
#' # Ewma estimation of the semideviation. lambda = 0.9
#' semidevEstimation(rets, control = list(type = 'ewma', lambda = 0.9))
#' @export
semidevEstimation <- function(rets, control = list()) {
  
  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x N) matrix")
  }
  
  ctr <- .ctrSemidev(control)
  
  if (ctr$type[1] == "naive") {
    semiDev <- .naiveSemiDev(rets = rets)
  } else if (ctr$type[1] == "ewma") {
    semiDev <- .ewmaSemiDev(rets = rets, lambda = ctr$lambda)
  } else {
    stop("control$type is not well defined")
  }
  return(semiDev)
}

.ctrSemidev <- function(control = list()) {
  
  if (!is.list(control)) {
    stop("control must be a list")
  }
  if (length(control) == 0) {
    control <- list(type = "naive", lambda = 0.94)
  }
  nam <- names(control)
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- c("naive", "ewma")
  }
  if (!("lambda" %in% nam) || is.null(control$lambda)) {
    control$lambda <- 0.94
  }
  return(control)
}

.naiveSemiDev <- function(rets) {
  ## Compute the naive semi-deviation INPUTs rets : matrix (T x N) returns
  ## OUTPUTs semiDev : vector (N x 1) semi-deviation
  semiDev <- .ewmaSemiDev(rets, lambda = 1)
  return(semiDev)
}

.ewmaSemiDev <- function(rets, lambda) {
  ## Compute the ewma semi-deviation INPUTs rets : matrix (T x N) returns
  ## OUTPUTs semiDev : vector (N x 1) semi-deviation
  t <- dim(rets)[1]
  n <- dim(rets)[2]
  semiDev <- vector("double", n)
  mu <- colMeans(rets)
  w <- lambda^(t:1)
  for (j in 1:n) {
    retsj <- rets[, j]
    muj <- mu[j]
    idx <- retsj < muj
    wj <- w[idx]/sum(w[idx])
    semiDev[j] <- sqrt(sum(wj * (retsj[idx] - muj)^2))
  }
  return(semiDev)
}