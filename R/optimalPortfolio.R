#' @name optimalPortfolio
#' @title Optimal portfolio
#' @description Function wich computes the optimal portfolio's weights.
#' @details The argument \code{control} is a list that can supply any of the following
#' components: 
#' \itemize{
#' \item \code{type} method used to compute the
#' optimal portfolio, among \code{'mv'}, \code{'minvol'}, \code{'invvol'},
#' \code{'erc'}, \code{'maxdiv'} and \code{'riskeff'} where: 
#' 
#' \code{'mv'} is used to compute the weights of the mean-variance portfolio. The weights are
#' computed following this equation: \deqn{w = \frac{1}{\gamma} \Sigma^{-1}
#' \mu }{ w = 1 / \gamma \Sigma^{-1} \mu}. 
#' 
#' \code{'minvol'} is used to compute the weights of the minimum variance portfolio.  
#' 
#' \code{'erc'} is used to compute the weights of the equal-risk-contribution portfolio. For a 
#' portfolio \eqn{w}, the percentage volatility risk contribution of the i-th
#' asset in the portfolio is given by: 
#' \deqn{\% RC_i = \frac{ w_i {[\Sigma w]}_i}{w' \Sigma w} }{ RC_i = w_i[\Sigma w]_i / (w' \Sigma w)}. 
#' Then we compute the optimal portfolio by solving the following optimization problem:
#' \deqn{w = argmin \left\{ \sum_{i=1}^N (\% RC_i - \frac{1}{N})^2 \right\}
#' }{ argmin { \sum_{i=1}^{N} (RC_i - 1/N)^2} }.
#' 
#' \code{'maxdiv'} is used to compute the weights of the maximum diversification portfolio where:
#' \deqn{DR(w) = \frac{ w' \sigma}{\sqrt{w' \Sigma w} } \geq 1 }{ DR(w) = (w'
#' \sigma)/(\sqrt(w' \Sigma w)) \ge 1} is used in the optimization problem.
#' 
#' \code{'riskeff'} is used to compute the weights of the risk-efficient
#' portfolio: \deqn{w = {argmax}\left\{ \frac{w' J \xi}{ \sqrt{w' \Sigma w}
#' }\right\} }{w = argmax ( w'J \xi)\sqrt(w'\Sigma w)} where \eqn{J} is a
#' \eqn{(N \times 10)}{(N x 10)} matrix of zeros whose \eqn{(i,j)}-th element
#' is one if the semi-deviation of stock \eqn{i} belongs to decile
#' \eqn{j},\eqn{\xi = (\xi_1,\ldots,\xi_{10})'}. 
#' 
#' \code{'invvol'} is the inverse volatility portfolio.
#' 
#' Default: \code{type = 'mv'}.
#' 
#' \item \code{constraint} constraint used for the optimization, among
#' \code{'none'}, \code{'lo'} and \code{'gross'} where: \code{'none'} is used to 
#' compute the unconstraint portfolio, \code{'lo'} is the long-only constraint and 
#' \code{'gross'} is the gross constraint. Default: \code{constraint = 'none'}. 
#' 
#' \item \code{gross.c} gross exposure constraint. Default: \code{gross.c = 1.6}. 
#' 
#' \item \code{gamma} Risk aversion parameter. Default: \code{gamma = 0.89}.
#' }
#' 
#' @param Sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param mu a \eqn{(N \times 1)}{(N x 1)} vector of expected returns. Default:
#' \code{mu = NULL}.
#' @param semiDev a vector \eqn{(N \times 1)}{(N x 1)} of semideviation.
#' Default: \code{semiDev = NULL}.
#' @param control control parameters (see *Details*).
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of optimal portfolio weights.
#' @note The long-only and the gross constraints are implemented for \code{'mv'}, \code{'minvol'} and 
#' \code{'maxdiv'} portfolios.
#' @author David Ardia, Kris Boudt and Jean-Philippe Gagnon Fleury.
#' @references Amenc, N., Goltz, F., Martellini, L., Retowsky, P. (2011).
#' Efficient indexation: An alternatice to cap-weightes indices.  \emph{Journal
#' of Investment Management} \bold{9}(4), pp.1--23.
#' 
#' Ardia, D., Boudt, K. (2015).  Implied expected returns and the choice of a
#' mean-variance efficient portfolio proxy.  \emph{Journal of Portfolio
#' Management} \bold{41} (4), pp.68--81.
#' 
#' Ardia, D., Bolliger, G., Boudt, K., Gagnon-Fleury, J.-P. (2016).  The Impact
#' of Covariance Misspecification in Risk-Based Portfolios.  \emph{Working
#' paper}. URL
#' \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2650644}
#' 
#' Choueifaty, Y., Coignard, Y., (2008).  Toward maximum diversification.
#' \emph{Journal of Portfolio Management} \bold{35} (1), pp.40--51.
#' 
#' Choueifaty, Y., Froidure, T., Reynier, J., (2011).  Properties of the most
#' diversified portfolio.  \emph{Working paper}.
#' 
#' Das, S., Markowitz, H., Scheid, J., Statman, M., (2010).  Portfolio
#' optimization with mental accounts.  \emph{Journal of Financial and
#' Quantitative Analysis} \bold{45}, pp.311--334.
#' 
#' DeMiguel, V., Garlappi, L., Uppal, R., (2009).  Optimal versus naive
#' diversification: How inefficient is the 1/n portfolio strategy.  \emph{The
#' Review of Financial Studies} \bold{22}(5), pp.1915--1953.
#' 
#' Fan, J., Zhang, J., Yu, K., March (2009).  Asset allocation and risk
#' assessment with gross exposure constraints for vast portfolios.
#' \emph{Working paper}.
#' 
#' Maillard, S., Roncalli, T., Teiletche, J., (2010).  The properties of
#' equally weighted risk contribution portfolios.  \emph{Journal of Portfolio
#' Management} \bold{36}(4), pp.60--70.
#' 
#' Martellini, L., (2008).  Towards the design of better equity benchmarks.
#' \emph{Journal of Portfolio Management} \bold{34}, Summer,pp.34--41.
#' @keywords optimize
#' @examples
#' # For the examples, we simply generate a 100 x 25 random matrix.
#' set.seed(3214)
#' T = 100
#' N = 25
#' rets = matrix(rnorm(T * N), nrow = T, ncol = N)
#' 
#' mu = meanEstimation(rets)
#' Sigma = covEstimation(rets)
#' semiDev = semidevEstimation(rets)
#' 
#' # Mean-variance portfolio without constraint and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma)
#' 
#' # Mean-variance portfolio without constraint and gamma = 1
#' optimalPortfolio(mu = mu, Sigma = Sigma, control = list(gamma = 1))
#' 
#' # Mean-variance portfolio with without and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv'))
#' 
#' # Mean-variance portfolio without constraint and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv', 
#'     constraint = 'none'))
#' 
#' # Mean-variance portfolio with the long-only constraint and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv', 
#'     constraint = 'lo'))
#' 
#' # Mean-variance portfolio with the gross constraint, gross constraint parameter = 1.6 and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv', 
#'     constraint = 'gross'))
#' 
#' # Mean-variance portfolio with the gross constraint, gross constraint parameter = 1.2 and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, control = list(type = 'mv', 
#'     constraint = 'gross', gross.c = 1.2))
#' 
#' # Minimum volatility portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol'))
#' 
#' # Minimum volatility portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol', constraint = 'none'))
#' 
#' # Minimim volatility portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol', constraint = 'lo'))
#' 
#' # Minimum volatility portfolio with the gross constraint and the gross constraint parameter = 1.6
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol', constraint = 'gross'))
#' 
#' # Minimum volatility portfolio with the gross constraint and the gross parameter = 1.2
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'minvol', constraint = 'gross', 
#'     gross.c = 1.2))
#' 
#' # Equal-risk-contribution portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'erc'))
#' 
#' # Equal-risk-contribution portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'erc', constraint = 'none'))
#' 
#' # Equal-risk-contribution portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'erc', constraint = 'lo'))
#' 
#' # Maximum diversification portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'maxdiv'))
#' 
#' # Maximum diversification portoflio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'maxdiv', constraint = 'lo'))
#' 
#' # Risk-efficient portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, semiDev = semiDev, control = list(type = 'riskeff'))
#' 
#' # Risk-efficient portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, semiDev = semiDev, control = list(type = 'riskeff', 
#'     constraint = 'lo'))
#'     
#' # Inverse volatility portfolio
#' optimalPortfolio(Sigma = Sigma, control = list(type = 'invvol'))
#' @export
optimalPortfolio <- function(Sigma, mu = NULL, semiDev = NULL, control = list()) {
  
  if (missing(Sigma)) {
    stop("A covariance matrix (Sigma) is required")
  }
  if (!is.matrix(Sigma)) {
    stop("Sigma must be a matrix")
  }
  if (!isSymmetric(Sigma)) {
    stop("Sigma must be a symmetric matrix")
  }
  
  ctr <- .ctrPortfolio(control)
  
  if (ctr$type[1] == "mv") {
    w <- .mvPortfolio(mu = mu, Sigma = Sigma, control = control)
  } else if (ctr$type[1] == "minvol") {
    w <- .minvolPortfolio(Sigma = Sigma, control = control)
  } else if (ctr$type[1] == "erc") {
    w <- .ercPortfolio(Sigma = Sigma, control = control)
  } else if (ctr$type[1] == "maxdiv") {
    w <- .maxdivPortfolio(Sigma = Sigma, control = control)
  } else if (ctr$type[1] == "riskeff") {
    w <- .riskeffPortfolio(Sigma = Sigma, semiDev = semiDev, control = control)
  } else if (ctr$type[1] == "invvol") {
    w <- .invvolPortfolio(Sigma = Sigma, control = control)
  } else {
    stop("control$type is not well defined")
  }
  return(w)
}

.ctrPortfolio <- function(control = list()) {
  ## Function used to control the list input INPUTs control : a control
  ## list The argument control is a list that can supply any of the
  ## following components type : 'mv', 'minvol', 'erc', 'maxdiv',
  ## riskeff' constraint : 'none', 'lo', 'gross' gross.c : default = 1.6
  ## gamma : default = 0.89 OUTPUTs control : list
  if (!is.list(control)) {
    stop("control must be a list")
  }
  if (length(control) == 0) {
    control <- list(type = "mv", constraint = "none", gross.c = 1.6, 
                    LB = NULL, UB = NULL, w0 = NULL, gamma = 0.89)
  }
  nam <- names(control)
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- c("mv", "minvol", "erc", "maxdiv", "riskeff", "invvol")
  }
  if (!("constraint" %in% nam) || is.null(control$constraint)) {
    control$constraint <- c("none", "lo", "gross")
  }
  if (!("gross.c" %in% nam) || is.null(control$gross.c)) {
    control$gross.c <- 1.6
  }
  if (!("LB" %in% nam) || is.null(control$LB)) {
    control$LB <- NULL
  }
  if (!("UB" %in% nam) || is.null(control$UB)) {
    control$UB <- NULL
  }
  if (!("w0" %in% nam) || is.null(control$w0)) {
    control$w0 <- NULL
  }
  if (!("gamma" %in% nam) || is.null(control$gamma)) {
    control$gamma <- c(0.8773, 2.7063, 3.795)
  }
  return(control)
}

.mvPortfolio <- function(mu, Sigma, control = list()) {
  ## Compute the weight of the mean-variance portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  ctr <- .ctrPortfolio(control)
  
  if (is.null(mu)) {
    stop("A vector of mean (mu) is required to compute the mean-variance portfolio")
  }
  if (ctr$constraint[1] == "none") {
    invSigmamu <- solve(Sigma, mu)
    w <- (1/ctr$gamma[1]) * invSigmamu/sum(invSigmamu)  
  } else if (ctr$constraint[1] == "lo") {
    n <- dim(Sigma)[1]
    Dmat <- ctr$gamma[1] * Sigma
    Amat <- cbind(rep(1, n), diag(n))
    bvec <- c(1, rep(0, n))
    w <- quadprog::solve.QP(Dmat = Dmat, dvec = mu, Amat = Amat, bvec = bvec, 
                            meq = 1)$solution
  } else if (ctr$constraint[1] == "gross") {
    .meanvar <- function(w) {
      Sigmaw <- crossprod(Sigma, w)
      opt <- -as.numeric(crossprod(mu, w)) + 0.5 * ctr$gamma[1] * 
        as.numeric(crossprod(w, Sigmaw))
      return(opt)
    }
    
    n <- dim(Sigma)[1]
    w0 <- ctr$w0
    if (is.null(w0)) {
      w0 <- rep(1/n, n)
    }
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
    w <- nloptr::slsqp(x0 = w0, fn = .meanvar, hin = ..grossContraint, 
                       heq = .eqConstraint, nl.info = FALSE, 
                       control = list(xtol_rel = 1e-18, check_derivatives = FALSE))$par
  } else {
    stop("control$constraint not well defined")
  }
  return(w)
}

.minvolPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the minimum volatility portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  ctr <- .ctrPortfolio(control)
  n <- dim(Sigma)[1]
  
  if (ctr$constraint[1] == "none") {
    tmp <- solve(Sigma, rep(1, n))
    w <- tmp/sum(tmp)
  } else if (ctr$constraint[1] == "lo") {
    dvec <- rep(0, n)
    Amat <- cbind(rep(1, n), diag(n))
    bvec <- c(1, rep(0, n))
    w <- quadprog::solve.QP(Dmat = Sigma, dvec = dvec, Amat = Amat, 
                            bvec = bvec, meq = 1)$solution
  } else if (ctr$constraint[1] == "gross") {
    .minvol <- function(w) {
      Sigmaw <- crossprod(Sigma, w)
      v <- as.numeric(crossprod(w, Sigmaw))
      return(v)
    }
    
    n <- dim(Sigma)[1]
    w0 <- ctr$w0
    if (is.null(w0)) {
      w0 <- rep(1/n, n)
    }
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
    w <- nloptr::slsqp(x0 = w0, fn = .minvol, hin = ..grossContraint, 
                       heq = .eqConstraint, nl.info = FALSE, 
                       control = list(xtol_rel = 1e-18, check_derivatives = FALSE))$par
  } else {
    stop("control$constraint not well defined")
  }
  return(w)
}

.ercPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the equal-risk-contribution portfolio INPUTs
  ## Sigma : matrix (N x N) covariance matrix control : list of control
  ## parameters OUTPUTs w : vector (N x 1) weight
  ctr <- .ctrPortfolio(control)
  
  n <- dim(Sigma)[2]
  w0 <- ctr$w0
  if (is.null(w0)) {
    # w0 = rep(1/n, n)
    w0 <- 1/sqrt(diag(Sigma))
    w0 <- w0/sum(w0)
  }
  
  .pRC <- function(w) {
    Sigmaw <- crossprod(Sigma, w)
    pRC <- (w * Sigmaw)/as.numeric(crossprod(w, Sigmaw))
    d <- sum((pRC - 1/n)^2)
    return(d)
  }
  
  # DA currently not implemented for long-only and gross constraints
  if (ctr$constraint[1] == "none") {
    w <- nloptr::slsqp(x0 = w0, fn = .pRC, heq = .eqConstraint, lower = rep(0, n), nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE))$par
  } else if (ctr$constraint[1] == "lo") {
    w <- nloptr::slsqp(x0 = w0, fn = .pRC, heq = .eqConstraint, lower = rep(0, n), nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else if (ctr$constraint[1] == "gross") {
    w <- nloptr::slsqp(x0 = w0, fn = .pRC, heq = .eqConstraint, lower = rep(0, n), nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else {
    stop("control$constraint not well defined")
  }
  return(w)
}

.maxdivPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of weight of the maximum diversification portfolio
  ## INPUTs Sigma : matrix (N x N) covariance matrix control : list of
  ## control parameters OUTPUTs w : vector (N x 1) weight
  ctr <- .ctrPortfolio(control)
  
  n <- dim(Sigma)[2]
  w0 <- ctr$w0
  if (is.null(w0)) {
    w0 <- rep(1/n, n)
  }
  
  .divRatio <- function(w) {
    sig <- sqrt(diag(Sigma))
    Sigmaw <- crossprod(Sigma, w)
    divRatio <- as.numeric(-crossprod(w, sig)/sqrt(crossprod(w, Sigmaw)))
    return(divRatio)
  }
  if (ctr$constraint[1] == "none") {
    w <- nloptr::slsqp(x0 = w0, fn = .divRatio, heq = .eqConstraint, 
                       nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else if (ctr$constraint[1] == "lo") {
    w <- nloptr::slsqp(x0 = w0, fn = .divRatio, lower = rep(0, n), 
                       heq = .eqConstraint, nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else if (ctr$constraint[1] == "gross") {
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
    w <- nloptr::slsqp(x0 = w0, fn = .divRatio, hin = ..grossContraint, 
                       heq = .eqConstraint, nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else {
    stop("control$constraint not well defined")
  }
  return(w)
}

.invvolPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the inverse-volatility portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  
  sig <- sqrt(diag(Sigma))
  w <- 1/sig
  w <- w/sum(w)
  return(w)
}

.riskeffPortfolio <- function(Sigma, semiDev, control = list()) {
  ## Compute the weight of the risk-efficient portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  ctr <- .ctrPortfolio(control)
  
  if (is.null(semiDev)) {
    stop("A vector of semideviation (semiDev) is require to compute the risk-efficient portfolio")
  }
  
  n <- dim(Sigma)[2]
  pct <- c(0, quantile(semiDev, probs = seq(0.1, 1, 0.1)))
  epsilon <- vector("double", n)
  J <- matrix(rep(0, n^2), ncol = n)
  
  for (i in 2:11) {
    pos <- semiDev > pct[i - 1] & semiDev <= pct[i]
    J[pos, i - 1] <- 1
    epsilon[i - 1] <- median(semiDev[pos])
  }
  Jepsilon <- crossprod(t(J), epsilon)
  # DA Additional constraints used to stabilize optimization
  LB <- (1/(2 * n)) * rep(1, n)
  UB <- (2/n) * rep(1, n)
  
  w0 <- ctr$w0
  if (is.null(w0)) {
    w0 <- (UB - LB)
    w0 <- w0/sum(w0)
  }
  
  .distRiskEff <- function(w) {
    Sigmaw <- crossprod(Sigma, w)
    d <- as.numeric(-crossprod(w, Jepsilon)/sqrt(crossprod(w, Sigmaw)))
    return(d)
  }
  
  .eqConstraintReff <- function(w) {
    return(sum(w))
  }
  
  if (ctr$constraint[1] == "none") {
    w <- nloptr::slsqp(x0 = w0, fn = .distRiskEff, heq = .eqConstraint, 
                       lower = LB, upper = UB, nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else if (ctr$constraint[1] == "lo") {
    w <- nloptr::slsqp(x0 = w0, fn = .distRiskEff, heq = .eqConstraint, 
                       lower = LB, upper = UB, nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else if (ctr$constraint[1] == "gross") {
    w <- nloptr::slsqp(x0 = w0, fn = .distRiskEff, heq = .eqConstraint, 
                       lower = LB, upper = UB, nl.info = FALSE, 
                       control = list(xtol_rel = 1e-08, check_derivatives = FALSE, maxeval = 2000))$par
  } else {
    stop("control$constraint not well defined")
  }
  return(w)
}

## Constraints used by the optimizers
.eqConstraint <- function(w) {
  return(sum(w) - 1)
}

# DA here 1.6 is hard-coded, this should be changed
.grossConstraint <- function(w, gross.c) {
  return(gross.c - norm(as.matrix(w), type = "1"))
}
