#' @name optimalPortfolio
#' @aliases optimalPortfolio
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
#' \mu }{w = 1 / \gamma \Sigma^{-1} \mu}. 
#' 
#' \code{'minvol'} is used to compute the weights of the minimum variance portfolio.  
#' 
#' \code{'invvol'} is the inverse volatility portfolio.
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
#' Default: \code{type = 'mv'}.
#' 
#' These portfolios are summarized in Ardia and Boudt (2015) and Ardia et al. (2016). Below we list the various references.
#' 
#' \item \code{constraint} constraint used for the optimization, among
#' \code{'none'}, \code{'lo'}, \code{'gross'} and \code{'user'}, where: \code{'none'} is used to 
#' compute the unconstraint portfolio, \code{'lo'} is the long-only constraints (non-negative weighted),  
#' \code{'gross'} is the gross exposure constraint, and \code{'user'} is the set of user constraints (typically
#' lower and upper boundaries. Default: \code{constraint = 'none'}. Note that the 
#' summability constraint is always imposed.
#' 
#' \item \code{LB} lower boundary for the weights. Default: \code{LB = NULL}. 
#' 
#' \item \code{UB} lower boundary for the weights. Default: \code{UB = NULL}. 
#' 
#' \item \code{w0} starting value for the optimizer. Default: \code{w0 = NULL} takes the 
#' equally-weighted portfolio as a starting value. When \code{LB} and \code{UB} are provided, it is set to 
#' mid-point of the bounds.
#' 
#' \item \code{gross.c} gross exposure constraint. Default: \code{gross.c = 1.6}. 
#' 
#' \item \code{gamma} risk aversion parameter. Default: \code{gamma = 0.89}.
#' 
#' \item \code{ctr.slsqp} list with control parameters for slsqp function.
#' }
#' 
#' @param Sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param mu a \eqn{(N \times 1)}{(N x 1)} vector of expected returns. Default:
#' \code{mu = NULL}.
#' @param semiDev a vector \eqn{(N \times 1)}{(N x 1)} of semideviations.
#' Default: \code{semiDev = NULL}.
#' @param control control parameters (see *Details*).
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of optimal portfolio weights.
#' @author David Ardia, Kris Boudt and Jean-Philippe Gagnon Fleury.
#' @references 
#' Amenc, N., Goltz, F., Martellini, L., Retowsky, P. (2011).
#' Efficient indexation: An alternatice to cap-weightes indices.  
#' \emph{Journal of Investment Management} \bold{9}(4), pp.1-23.
#' 
#' Ardia, D., Boudt, K. (2015). 
#' Implied expected returns and the choice of a mean-variance efficient portfolio proxy.
#' \emph{Journal of Portfolio Management} \bold{41}(4), pp.66-81. 
#' \doi{10.3905/jpm.2015.41.4.068}
#' 
#' Ardia, D., Bolliger, G., Boudt, K., Gagnon-Fleury, J.-P. (2016).  
#' \emph{The Impact of covariance misspecification in risk-based portfolios}.  
#' Working paper. 
#' \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2650644}
#' 
#' Choueifaty, Y., Coignard, Y. (2008).  
#' Toward maximum diversification.
#' \emph{Journal of Portfolio Management} \bold{35}(1), pp.40-51. 
#' \doi{10.3905/JPM.2008.35.1.40}
#' 
#' Choueifaty, Y., Froidure, T., Reynier, J. (2013).  
#' Properties of the most diversified portfolio.  
#' \emph{Journal of Investment Strategies} \bold{2}(2), pp.49-70. 
#' \doi{10.21314/JOIS.2013.033}
#' 
#' Das, S., Markowitz, H., Scheid, J., Statman, M. (2010).  
#' Portfolio optimization with mental accounts.  
#' \emph{Journal of Financial and Quantitative Analysis} \bold{45}(2), pp.311-334. 
#' \doi{10.1017/S0022109010000141}
#' 
#' DeMiguel, V., Garlappi, L., Uppal, R. (2009).  
#' Optimal versus naive diversification: How inefficient is the 1/n portfolio strategy.  
#' \emph{Review of Financial Studies} \bold{22}(5), pp.1915-1953. 
#' \doi{10.1093/rfs/hhm075}
#' 
#' Fan, J., Zhang, J., Yu, K. (2009).  
#' \emph{Asset allocation and risk assessment with gross exposure constraints for vast portfolios}.
#' Working paper. 
#' \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1307423}
#' 
#' Maillard, S., Roncalli, T., Teiletche, J. (2010).  
#' The properties of equally weighted risk contribution portfolios.  
#' \emph{Journal of Portfolio Management} \bold{36}(4), pp.60-70. 
#' \doi{10.3905/jpm.2010.36.4.060}
#' 
#' Martellini, L. (2008).  
#' Towards the design of better equity benchmarks.
#' \emph{Journal of Portfolio Management} \bold{34}(4), Summer,pp.34-41. 
#' \doi{10.3905/jpm.2008.709978}
#' @keywords optimize
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets = Industry_10
#' 
#' # Mean estimation
#' mu = meanEstimation(rets)
#' 
#' # Covariance estimation
#' Sigma = covEstimation(rets)
#' 
#' # Semi-deviation estimation
#' semiDev = semidevEstimation(rets)
#' 
#' # Mean-variance portfolio without constraint and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma)
#' 
#' # Mean-variance portfolio without constraint and gamma = 1
#' optimalPortfolio(mu = mu, Sigma = Sigma, 
#'   control = list(gamma = 1))
#' 
#' # Mean-variance portfolio without constraint and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, 
#'   control = list(type = 'mv'))
#' 
#' # Mean-variance portfolio without constraint and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, 
#'   control = list(type = 'mv', constraint = 'none'))
#' 
#' # Mean-variance portfolio with the long-only constraint and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, 
#'   control = list(type = 'mv', constraint = 'lo'))
#' 
#' # Mean-variance portfolio with LB and UB constraints
#' optimalPortfolio(mu = mu, Sigma = Sigma, 
#'   control = list(type = 'mv', constraint = 'user', LB = rep(0.02, 10), UB = rep(0.8, 10)))
#' 
#' # Mean-variance portfolio with the gross constraint, 
#' # gross constraint parameter = 1.6 and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, 
#'   control = list(type = 'mv', constraint = 'gross'))
#' 
#' # Mean-variance portfolio with the gross constraint, 
#' # gross constraint parameter = 1.2 and gamma = 0.89
#' optimalPortfolio(mu = mu, Sigma = Sigma, 
#'   control = list(type = 'mv', constraint = 'gross', gross.c = 1.2))
#' 
#' # Minimum volatility portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'minvol'))
#' 
#' # Minimum volatility portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'minvol', constraint = 'none'))
#' 
#' # Minimim volatility portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'minvol', constraint = 'lo'))
#'   
#' # Minimim volatility portfolio with LB and UB constraints
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'minvol', constraint = 'user', LB = rep(0.02, 10), UB = rep(0.8, 10)))
#' 
#' # Minimum volatility portfolio with the gross constraint 
#' # and the gross constraint parameter = 1.6
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'minvol', constraint = 'gross'))
#' 
#' # Minimum volatility portfolio with the gross constraint 
#' # and the gross parameter = 1.2
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'minvol', constraint = 'gross', gross.c = 1.2))
#'     
#' # Inverse volatility portfolio
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'invvol'))
#' 
#' # Equal-risk-contribution portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'erc', constraint = 'lo'))
#'   
#' # Equal-risk-contribution portfolio with LB and UB constraints
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'erc', constraint = 'user', LB = rep(0.02, 10), UB = rep(0.8, 10)))
#' 
#' # Maximum diversification portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'maxdiv'))
#' 
#' # Maximum diversification portoflio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'maxdiv', constraint = 'lo'))
#'   
#' # Maximum diversification portoflio with LB and UB constraints
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'maxdiv', constraint = 'user', LB = rep(0.02, 10), UB = rep(0.8, 10)))
#' 
#' # Risk-efficient portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, semiDev = semiDev, 
#'   control = list(type = 'riskeff'))
#' 
#' # Risk-efficient portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, semiDev = semiDev, 
#'   control = list(type = 'riskeff', constraint = 'lo'))
#'   
#' # Risk-efficient portfolio with LB and UB constraints
#' optimalPortfolio(Sigma = Sigma, semiDev = semiDev, 
#'   control = list(type = 'riskeff', constraint = 'user', LB = rep(0.02, 10), UB = rep(0.8, 10)))
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
  
  n = dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)
  
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

.ctrPortfolio <- function(n, control = list()) {
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
                    LB = NULL, UB = NULL, w0 = NULL, gamma = 0.89, ctr.slsqp = NULL)
  }
  nam <- names(control)
  ## type
  type <- c("mv", "minvol", "erc", "maxdiv", "riskeff", "invvol")
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- type
  }
  if (!(control$type[1] %in% type)) {
    stop("'type' is not properly defined")
  }
  
  ## constraint
  constraint <- c("none", "lo", "gross", "user")
  if (!("constraint" %in% nam) || is.null(control$constraint)) {
    control$constraint <- constraint
  }
  if (!(control$constraint[1] %in% constraint)) {
    stop("'constraint' is not properly defined")
  }
  
  ## gross.c
  if (!("gross.c" %in% nam) || is.null(control$gross.c)) {
    control$gross.c <- 1.6
  }
  if (constraint[1] == "gross") {
    if (!is.numeric(control$gross.c)) {
      stop("'gross.c' is not properly defined")
    }
  }
  
  ## user LB and UB
  if (!("LB" %in% nam) || is.null(control$LB)) {
    control$LB <- NULL
  }
  if (!("UB" %in% nam) || is.null(control$UB)) {
    control$UB <- NULL
  }
  if (control$constraint[1] == "user") {
    if (is.null(control$LB) || (length(control$LB) != n)) {
      stop("'LB' is not properly defined")
    }
    if (is.null(control$UB) || (length(control$UB) != n)) {
      stop("'UB' is not properly defined")
    }
  }
  if (control$constraint[1] == "lo") {
    control$LB <- rep(0, n)
  }
  
  ## starting portfolio
  if (!("w0" %in% nam) || is.null(control$w0)) {
    control$w0 <- rep(1, n) / n
    if (!is.null(control$LB) && !is.null(control$UB)) {
      control$w0 = 0.5 * (control$LB + control$UB)
    }
  }
  if (length(control$w0) != n) {
    stop("'w0' is not properly defined")
  }
  
  # risk aversion parameter
  if (!("gamma" %in% nam) || is.null(control$gamma)) {
    control$gamma <- c(0.8773, 2.7063, 3.795)
  }
  if (!is.numeric(control$gamma)) {
    stop("'gamma' is not properly defined")
  }
  
  # optimization list
  if (!("ctr.slsqp" %in% nam) || is.null(control$ctr.slsqp)) {
    control$ctr.slsqp <- list(xtol_rel = 1e-18, check_derivatives = FALSE, maxeval = 2000)
  }
  if (!is.list(control$ctr.slsqp)) {
    stop("'ctr.slsqp' is not properly defined")
  }
  
  return(control)
}

.mvPortfolio <- function(mu, Sigma, control = list()) {
  ## Compute the weight of the mean-variance portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)
  
  if (is.null(mu)) {
    stop("A vector of mean (mu) is required to compute the mean-variance portfolio")
  }
  if (ctr$constraint[1] == "none") {
    invSigmamu <- solve(Sigma, mu)
    w <- (1/ctr$gamma[1]) * invSigmamu/sum(invSigmamu) 
  } else if (ctr$constraint[1] == "lo" || ctr$constraint[1] == "user") {
    Dmat <- ctr$gamma[1] * Sigma
    Amat <- cbind(rep(1, n), diag(n))
    bvec <- c(1, ctr$LB)
    if (ctr$constraint[1] == "user") {
      Amat <- cbind(Amat, -diag(n))
      bvec <- c(bvec, -ctr$UB)
    }
    w <- quadprog::solve.QP(Dmat = Dmat, dvec = mu, Amat = Amat, bvec = bvec, 
                            meq = 1)$solution
  } else if (ctr$constraint[1] == "gross") {
    .meanvar <- function(w) {
      Sigmaw <- crossprod(Sigma, w)
      opt <- -as.numeric(crossprod(mu, w)) + 0.5 * ctr$gamma[1] * 
        as.numeric(crossprod(w, Sigmaw))
      return(opt)
    }
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
    w <- nloptr::slsqp(x0 = ctr$w0, fn = .meanvar, 
                       hin = ..grossContraint, 
                       heq = .eqConstraint, 
                       lower = ctr$LB,
                       upper = ctr$UB,
                       nl.info = FALSE, control = ctr$ctr.slsqp)$par
  } else {
    # spotted in controls
  }
  w[w<=ctr$LB] <- ctr$LB[w<=ctr$LB]
  w[w>=ctr$UB] <- ctr$UB[w>=ctr$UB]
  w <- w / sum(w)
  return(w)
}

.minvolPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the minimum volatility portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)

  if (ctr$constraint[1] == "none") {
    tmp <- solve(Sigma, rep(1, n))
    w <- tmp/sum(tmp)
  } else if (ctr$constraint[1] == "lo" || ctr$constraint[1] == "user") {
    dvec <- rep(0, n)
    Amat <- cbind(rep(1, n), diag(n))
    bvec <- c(1, ctr$LB)
    if (ctr$constraint[1] == "user") {
      Amat <- cbind(Amat, -diag(n))
      bvec <- c(bvec, -ctr$UB)
    }
    w <- quadprog::solve.QP(Dmat = Sigma, dvec = dvec, Amat = Amat, 
                            bvec = bvec, meq = 1)$solution
  } else if (ctr$constraint[1] == "gross") {
    .minvol <- function(w) {
      Sigmaw <- crossprod(Sigma, w)
      v <- as.numeric(crossprod(w, Sigmaw))
      return(v)
    }
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
    w <- nloptr::slsqp(x0 = ctr$w0, fn = .minvol, 
                       hin = ..grossContraint, 
                       heq = .eqConstraint,
                       lower = ctr$LB,
                       upper = ctr$UB,
                       nl.info = FALSE,  control = ctr$ctr.slsqp)$par
  } else {
    # spotted in controls
  }
  w[w<=ctr$LB] <- ctr$LB[w<=ctr$LB]
  w[w>=ctr$UB] <- ctr$UB[w>=ctr$UB]
  w <- w / sum(w)
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

.ercPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the equal-risk-contribution portfolio INPUTs
  ## Sigma : matrix (N x N) covariance matrix control : list of control
  ## parameters OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[2]
  ctr <- .ctrPortfolio(n, control)
  # DA user could use this instead of equal weighted
  # w0 <- 1/sqrt(diag(Sigma))
  # w0 <- w0/sum(w0)
  
  .pRC <- function(w) {
    Sigmaw <- crossprod(Sigma, w)
    pRC <- (w * Sigmaw)/as.numeric(crossprod(w, Sigmaw))
    d <- sum((pRC - 1/n)^2)
    return(d)
  }
  
  ..grossContraint = NULL
  if (ctr$constraint[1] == "gross") {
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
  }
  
  w <- nloptr::slsqp(x0 = ctr$w0, fn = .pRC, 
                     hin = ..grossContraint,
                     heq = .eqConstraint, 
                     lower = ctr$LB, 
                     upper = ctr$UB, 
                     nl.info = FALSE, control = ctr$ctr.slsqp)$par
  
  w[w<=ctr$LB] <- ctr$LB[w<=ctr$LB]
  w[w>=ctr$UB] <- ctr$UB[w>=ctr$UB]
  w <- w / sum(w)
  return(w)
}

.maxdivPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of weight of the maximum diversification portfolio
  ## INPUTs Sigma : matrix (N x N) covariance matrix control : list of
  ## control parameters OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[2]
  ctr <- .ctrPortfolio(n, control)

  .divRatio <- function(w) {
    sig <- sqrt(diag(Sigma))
    Sigmaw <- crossprod(Sigma, w)
    divRatio <- as.numeric(-crossprod(w, sig)/sqrt(crossprod(w, Sigmaw)))
    return(divRatio)
  }
  
  ..grossContraint = NULL
  if (ctr$constraint[1] == "gross") {
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
  }
  
  w <- nloptr::slsqp(x0 = ctr$w0, fn = .divRatio, 
                     hin = ..grossContraint,
                     heq = .eqConstraint, 
                     lower = ctr$LB, 
                     upper = ctr$UB, 
                     nl.info = FALSE, control = ctr$ctr.slsqp)$par
  
  w[w<=ctr$LB] <- ctr$LB[w<=ctr$LB]
  w[w>=ctr$UB] <- ctr$UB[w>=ctr$UB]
  w <- w / sum(w)
  return(w)
}

.riskeffPortfolio <- function(Sigma, semiDev, control = list()) {
  ## Compute the weight of the risk-efficient portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[2]
  ctr <- .ctrPortfolio(n, control)
  
  if (is.null(semiDev)) {
    stop("A vector of semideviations (semiDev) is require to compute the risk-efficient portfolio")
  }
  
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
  LB <- ctr$LB
  if (is.null(LB)) {
    LB <- (1/(2 * n)) * rep(1, n)
  }
  UB <- ctr$UB
  if (is.null(UB)) {
    UB <- (2/n) * rep(1, n)
  }

  # DA overwrite w0 for better starting values
  w0 <- (UB - LB)
  w0 <- w0/sum(w0)
  
  .distRiskEff <- function(w) {
    Sigmaw <- crossprod(Sigma, w)
    d <- as.numeric(-crossprod(w, Jepsilon)/sqrt(crossprod(w, Sigmaw)))
    return(d)
  }
  
  ..grossContraint = NULL
  if (ctr$constraint[1] == "gross") {
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
  }
  
  w <- nloptr::slsqp(x0 = ctr$w0, fn = .distRiskEff, 
                     hin = ..grossContraint,
                     heq = .eqConstraint, 
                     lower = LB, 
                     upper = UB, 
                     nl.info = FALSE, control = ctr$ctr.slsqp)$par
  
  w[w<=ctr$LB] <- ctr$LB[w<=ctr$LB]
  w[w>=ctr$UB] <- ctr$UB[w>=ctr$UB]
  w <- w / sum(w)
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
