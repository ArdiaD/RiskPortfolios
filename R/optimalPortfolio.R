#' @name optimalPortfolio
#' @aliases optimalPortfolio
#' @title Optimal portfolio
#' @description Function which computes the optimal portfolio's weights.
#' @details The argument \code{control} is a list that can supply any of the following
#' components: 
#' \itemize{
#' \item \code{type} method used to compute the
#' optimal portfolio, among \code{'mv'}, \code{'minvol'}, \code{'invvol'},
#' \code{'erc'}, \code{'maxdiv'}, \code{'riskeff'} and \code{'maxdec'} where: 
#' 
#' Every supplied bound enters the optimization itself. Under the \code{'gross'}
#' constraint the mean-variance, minimum-variance and maximum-decorrelation
#' problems are recast as quadratic programs in \eqn{(w^+, w^-)}{(w+, w-)} with
#' \eqn{w = w^+ - w^-}{w = w+ - w-}, which turns \eqn{\|w\|_1 \le c}{||w||_1 <= c}
#' into a linear constraint and makes them exactly solvable.
#'
#' \code{'mv'} is used to compute the weights of the mean-variance portfolio, that
#' is the solution of \deqn{\max_w \left\{ w' \mu - \frac{\gamma}{2} w' \Sigma w \right\}
#' \quad s.t. \quad w'1 = 1}{max_w {w'mu - gamma/2 w'Sigma w} s.t. w'1 = 1}
#' and, in the unconstrained case, is given in two-fund form by
#' \deqn{w = \frac{\Sigma^{-1} 1}{1' \Sigma^{-1} 1} + \frac{1}{\gamma} \left( \Sigma^{-1} \mu -
#' \frac{1' \Sigma^{-1} \mu}{1' \Sigma^{-1} 1} \Sigma^{-1} 1 \right)}{w = Sigma^-1 1 / (1'Sigma^-1 1)
#' + 1/gamma (Sigma^-1 mu - (1'Sigma^-1 mu)/(1'Sigma^-1 1) Sigma^-1 1)}
#' that is, the minimum-variance portfolio plus \eqn{1/\gamma}{1/gamma} times a
#' self-financing speculative portfolio. Note that the summability constraint is
#' imposed here as it is for the other portfolios, so the risk aversion
#' \eqn{\gamma}{gamma} governs the tilt away from the minimum-variance portfolio
#' rather than the leverage.
#'
#' \code{'minvol'} is used to compute the weights of the minimum variance portfolio.
#'
#' \code{'invvol'} is the inverse volatility portfolio. It is computed in closed
#' form and is always long-only with weights summing to one, so the \code{'lo'}
#' and \code{'gross'} constraints hold automatically. Bounds do not enter its
#' definition: if \code{LB}/\code{UB} are supplied and bind, the weights are
#' projected onto \eqn{\{LB \le w \le UB, w'1 = 1\}}{{LB <= w <= UB, w'1 = 1}}
#' and a warning is issued.
#'
#' \code{'erc'} is used to compute the weights of the equal-risk-contribution portfolio. For a
#' portfolio \eqn{w}, the percentage volatility risk contribution of the i-th
#' asset in the portfolio is given by:
#' \deqn{\% RC_i = \frac{ w_i {[\Sigma w]}_i}{w' \Sigma w} }{ RC_i = w_i[\Sigma w]_i / (w' \Sigma w)}.
#' Then we compute the optimal portfolio by solving the following optimization problem:
#' \deqn{w = argmin \left\{ \sum_{i=1}^N (\% RC_i - \frac{1}{N})^2 \right\}
#' }{ argmin { \sum_{i=1}^{N} (RC_i - 1/N)^2} }.
#' Without bounds, or under the long-only constraint, the solution is obtained
#' exactly by solving \eqn{y_i [\Sigma y]_i = 1/N}{y[i] [Sigma y][i] = 1/N} for
#' \eqn{y > 0}{y > 0} and rescaling to \eqn{w = y / (y'1)}{w = y / (y'1)}, which
#' attains a zero objective. See Spinu (2013) and Chaves et al. (2012). This is
#' necessary because the objective above is not convex and admits local minima on
#' the boundary at which whole groups of assets receive a zero weight; a
#' general-purpose optimizer started at the equally-weighted portfolio converges
#' to them whenever the correlation matrix has sizeable negative entries. Under
#' \code{'user'} or \code{'gross'} constraints the problem is solved numerically
#' from several starting values.
#'
#' \code{'maxdiv'} is used to compute the weights of the maximum diversification portfolio where:
#' \deqn{DR(w) = \frac{ w' \sigma}{\sqrt{w' \Sigma w} } \geq 1 }{ DR(w) = (w'
#' \sigma)/(\sqrt(w' \Sigma w)) \ge 1} is used in the optimization problem.
#' 
#' \code{'riskeff'} is used to compute the weights of the risk-efficient
#' portfolio: \deqn{w = {argmax}\left\{ \frac{w' J \xi}{ \sqrt{w' \Sigma w}
#' }\right\} }{w = argmax (w'J xi) / sqrt(w'Sigma w)} where \eqn{J} is a
#' \eqn{(N \times 10)}{(N x 10)} matrix of zeros whose \eqn{(i,j)}-th element
#' is one if the semi-deviation of stock \eqn{i} belongs to decile
#' \eqn{j},\eqn{\xi = (\xi_1,\ldots,\xi_{10})'}.
#'
#' Note that this portfolio carries additional bounds that stabilize the
#' optimization: whenever \code{LB} (resp. \code{UB}) is not supplied it is set
#' to \eqn{1/(2N)}{1/(2N)} (resp. \eqn{2/N}{2/N}), under every value of
#' \code{constraint}. The risk-efficient portfolio is therefore never the
#' summability-only problem, it is always long-only, and \code{gross.c} can
#' never bind on it since \eqn{\|w\|_1 = 1}{||w||_1 = 1}. Supply \code{LB} and
#' \code{UB} explicitly to override these defaults. 
#' 
#' \code{'maxdec'} is used to compute the weights of the maximum-decorrelation
#' portfolio: \deqn{w = {argmax}\left\{ 1 -  \sqrt{w' R w} \right\}
#' }{w = argmax {1- \sqrt(w' R w)}} where \eqn{R} is the correlation matrix.
#' 
#' Default: \code{type = 'mv'}.
#' 
#' These portfolios are summarized in Ardia and Boudt (2015) and Ardia et al. (2017). Below we list the various references.
#' 
#' \item \code{constraint} constraint used for the optimization, among
#' \code{'none'}, \code{'lo'}, \code{'gross'} and \code{'user'}, where: \code{'none'} is used to
#' compute the unconstraint portfolio, \code{'lo'} is the long-only constraints (non-negative weighted),
#' \code{'gross'} is the gross exposure constraint, and \code{'user'} is the set of user constraints (typically
#' lower and upper boundaries). Default: \code{constraint = 'none'}. Note that the
#' summability constraint is always imposed.
#'
#' \item \code{LB} lower boundary for the weights. Default: \code{LB = NULL}.
#'
#' \item \code{UB} upper boundary for the weights. Default: \code{UB = NULL}.
#'
#' \item \code{w0} starting value for the optimizer. Default: \code{w0 = NULL} takes the
#' equally-weighted portfolio as a starting value. When \code{LB} and \code{UB} are provided, it is set to
#' mid-point of the bounds. For the non-convex problems (\code{'erc'} under
#' \code{'user'}/\code{'gross'} constraints, \code{'maxdiv'} and \code{'riskeff'})
#' \code{w0} is one of several starting values that are tried, the best solution
#' being returned.
#'
#' \item \code{gross.c} gross exposure constraint, that is the bound imposed on
#' \eqn{\|w\|_1}{||w||_1}. Since the weights sum to one, it cannot be smaller
#' than one. Default: \code{gross.c = 1.6}.
#'
#' \item \code{gamma} risk aversion parameter, a single positive number.
#' Default: \code{gamma = 0.89}.
#'
#' \item \code{ctr.slsqp} list with control parameters for the
#' \code{\link[nloptr]{slsqp}} function. Default:
#' \code{list(xtol_rel = 1e-18, ftol_rel = 1e-12, check_derivatives = FALSE,
#' maxeval = 2000)}. \code{xtol_rel} is below double precision and so cannot be
#' met on its own; \code{ftol_rel} supplies an attainable stopping criterion, and
#' the optimizer warns if it stops on \code{maxeval} instead of converging. It is
#' used only by the non-convex problems -- \code{'erc'} under \code{'user'} or
#' \code{'gross'} constraints, \code{'maxdiv'} and \code{'riskeff'}. The
#' mean-variance, minimum-variance and maximum-decorrelation portfolios are
#' quadratic programs under every constraint, including \code{'gross'}, and are
#' solved exactly by \code{\link[quadprog]{solve.QP}}.
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
#' Efficient indexation: An alternative to cap-weighted indices.  
#' \emph{Journal of Investment Management} \bold{9}(4), pp.1-23.
#' 
#' Ardia, D., Boudt, K. (2015). 
#' Implied expected returns and the choice of a mean-variance efficient portfolio proxy.
#' \emph{Journal of Portfolio Management} \bold{41}(4), pp.66-81.
#' \doi{10.3905/jpm.2015.41.4.068} 
#' 
#' Ardia, D., Bolliger, G., Boudt, K., Gagnon-Fleury, J.-P. (2017).  
#' The Impact of covariance misspecification in risk-based portfolios.  
#' \emph{Annals of Operations Research} \bold{254}(1-2), pp.1-16. 
#' \doi{10.1007/s10479-017-2474-7}
#' 
#' Chaves, D., Hsu, J., Li, F., Shakernia, O. (2012).
#' Efficient algorithms for computing risk parity portfolio weights.
#' \emph{Journal of Investing} \bold{21}(3), pp.150-163.
#'
#' Choueifaty, Y., Coignard, Y. (2008).
#' Toward maximum diversification.
#' \emph{Journal of Portfolio Management} \bold{35}(1), pp.40-51.
#'
#' Choueifaty, Y., Froidure, T., Reynier, J. (2013).
#' Properties of the most diversified portfolio.
#' \emph{Journal of Investment Strategies} \bold{2}(2), pp.49-70.
#'
#' Das, S., Markowitz, H., Scheid, J., Statman, M. (2010).  
#' Portfolio optimization with mental accounts.  
#' \emph{Journal of Financial and Quantitative Analysis} \bold{45}(2), pp.311-334. 
#' 
#' DeMiguel, V., Garlappi, L., Uppal, R. (2009).  
#' Optimal versus naive diversification: How inefficient is the 1/n portfolio strategy.  
#' \emph{Review of Financial Studies} \bold{22}(5), pp.1915-1953. 
#' 
#' Fan, J., Zhang, J., Yu, K. (2012).  
#' Vast portfolio selection with gross-exposure constraints.
#' \emph{Journal of the American Statistical Association} \bold{107}(498), pp.592-606. 
#' 
#' Maillard, S., Roncalli, T., Teiletche, J. (2010).  
#' The properties of equally weighted risk contribution portfolios.  
#' \emph{Journal of Portfolio Management} \bold{36}(4), pp.60-70. 
#' 
#' Martellini, L. (2008).
#' Towards the design of better equity benchmarks.
#' \emph{Journal of Portfolio Management} \bold{34}(4), Summer,pp.34-41.
#'
#' Spinu, F. (2013).
#' An algorithm for computing risk parity weights.
#' \emph{SSRN working paper}.
#' \doi{10.2139/ssrn.2297383}
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
#' # Maximum diversification portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'maxdiv', constraint = 'lo'))
#'   
#' # Maximum diversification portfolio with LB and UB constraints
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
#'   
#' # Maximum decorrelation portfolio without constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'maxdec'))
#' 
#' # Maximum decorrelation portfolio with the long-only constraint
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'maxdec', constraint = 'lo'))
#'   
#' # Maximum decorrelation portfolio with LB and UB constraints
#' optimalPortfolio(Sigma = Sigma, 
#'   control = list(type = 'maxdec', constraint = 'user', LB = rep(0.02, 10), UB = rep(0.8, 10)))
#' @export
#' @importFrom stats cov2cor
optimalPortfolio <- function(Sigma, mu = NULL, semiDev = NULL, control = list()) {
  
  if (missing(Sigma)) {
    stop("A covariance matrix (Sigma) is required")
  }
  if (!is.matrix(Sigma)) {
    stop("Sigma must be a matrix")
  }
  if (!is.numeric(Sigma) || any(!is.finite(Sigma))) {
    stop("Sigma must be numeric and must not contain missing or infinite values")
  }
  if (!isSymmetric(Sigma)) {
    stop("Sigma must be a symmetric matrix")
  }
  if (any(diag(Sigma) < 0)) {
    stop("Sigma must have non-negative diagonal entries")
  }
  .checkPSD(Sigma)

  n = dim(Sigma)[1]
  if (!is.null(mu) && (length(mu) != n || !is.numeric(mu) || any(!is.finite(mu)))) {
    stop("mu must be a finite numeric vector of the same dimension as Sigma")
  }
  if (!is.null(semiDev) && (length(semiDev) != n || !is.numeric(semiDev) ||
                            any(!is.finite(semiDev)))) {
    stop("semiDev must be a finite numeric vector of the same dimension as Sigma")
  }
  if (!is.null(semiDev) && any(semiDev < 0)) {
    stop("semiDev must be non-negative")
  }
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
  } else if (ctr$type[1] == "maxdec") {
    w <- .maxdecPortfolio(Sigma = Sigma, control = control)
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
  nam <- names(control)
  ## type
  type <- c("mv", "minvol", "erc", "maxdiv", "riskeff", "invvol", "maxdec")
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
  if (control$constraint[1] == "gross") {
    if (!is.numeric(control$gross.c) || length(control$gross.c) != 1L ||
        !is.finite(control$gross.c)) {
      stop("'gross.c' is not properly defined")
    }
    if (control$gross.c < 1) {
      stop("'gross.c' must be at least one: the weights sum to one, so ",
           "||w||_1 >= 1 for any admissible portfolio")
    }
  }

  ## user LB and UB
  if (!("LB" %in% nam) || is.null(control$LB)) {
    control$LB <- NULL
  }
  if (!("UB" %in% nam) || is.null(control$UB)) {
    control$UB <- NULL
  }
  if (control$constraint[1] == "user" && (is.null(control$LB) || is.null(control$UB))) {
    stop("'LB' and 'UB' are required when constraint = 'user'")
  }
  ## the bounds are checked whenever they are supplied, not only under the
  ## 'user' constraint: they are passed to the optimizer under 'gross' as well,
  ## and bounds incompatible with 1'w = 1 make the problem infeasible
  if (!is.null(control$LB) && ((length(control$LB) != n) || !all(is.finite(control$LB)))) {
    stop("'LB' is not properly defined")
  }
  if (!is.null(control$UB) && ((length(control$UB) != n) || !all(is.finite(control$UB)))) {
    stop("'UB' is not properly defined")
  }
  if (!is.null(control$LB) && !is.null(control$UB) && any(control$LB > control$UB)) {
    stop("'LB' must not exceed 'UB'")
  }
  if ((!is.null(control$LB) && sum(control$LB) > 1) ||
      (!is.null(control$UB) && sum(control$UB) < 1)) {
    stop("'LB' and 'UB' are incompatible with the summability constraint")
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
  if (length(control$w0) != n || !all(is.finite(control$w0))) {
    stop("'w0' is not properly defined")
  }
  ## slsqp rejects a starting value outside the bounds, so w0 -- whether it is
  ## the equally-weighted default or was supplied by the caller -- is projected
  ## onto {LB <= w <= UB, w'1 = 1}. Without this, a perfectly feasible problem
  ## such as constraint = 'gross' with LB = c(0.9, 0, 0, 0) failed outright,
  ## because the default w0 = 1/n violates LB.
  if (!is.null(control$LB) || !is.null(control$UB)) {
    control$w0 <- .projectBox(control$w0, control$LB, control$UB)
  }
  ## the gross budget applies to the starting value too: slsqp handed a start
  ## that already violates it can stop there and return it unchanged
  if (control$constraint[1] == "gross" && sum(abs(control$w0)) > control$gross.c + 1e-08) {
    control$w0 <- .projectBox(rep(1, n)/n, control$LB, control$UB)
    if (sum(abs(control$w0)) > control$gross.c + 1e-08) {
      stop("'gross.c' is incompatible with 'LB' and 'UB'")
    }
  }

  # risk aversion parameter
  if (!("gamma" %in% nam) || is.null(control$gamma)) {
    control$gamma <- 0.89
  }
  if (!is.numeric(control$gamma) || length(control$gamma) != 1L ||
      !is.finite(control$gamma) || control$gamma <= 0) {
    stop("'gamma' must be a single positive number")
  }

  # optimization list
  ## xtol_rel alone is below double precision and can essentially never be met,
  ## so without ftol_rel the optimizer runs to maxeval on every gross-constrained
  ## problem and its convergence status carries no information
  if (!("ctr.slsqp" %in% nam) || is.null(control$ctr.slsqp)) {
    control$ctr.slsqp <- list(xtol_rel = 1e-18, ftol_rel = 1e-12,
                              check_derivatives = FALSE, maxeval = 2000)
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
  bounded <- !is.null(ctr$LB) || !is.null(ctr$UB)
  if (ctr$constraint[1] == "gross") {
    w <- .grossQP(ctr$gamma[1] * Sigma, mu, ctr, n)
  } else if (!bounded) {
    ## Solution of max_w {w'mu - gamma/2 w'Sigma w} s.t. w'1 = 1, in two-fund
    ## form: the minimum-variance portfolio plus 1/gamma times a self-financing
    ## speculative portfolio. Rescaling Sigma^-1 mu to sum to one instead would
    ## cancel gamma altogether, and would flip the sign of the whole portfolio
    ## whenever 1'Sigma^-1 mu < 0.
    invSigma1 <- solve(Sigma, rep(1, n))
    invSigmamu <- solve(Sigma, mu)
    a <- sum(invSigma1)
    b <- sum(invSigmamu)
    w <- invSigma1/a + (1/ctr$gamma[1]) * (invSigmamu - (b/a) * invSigma1)
  } else {
    w <- .boundedQP(ctr$gamma[1] * Sigma, mu, ctr, n)
  }
  w <- .finalizeWeights(w, ctr)
  .checkGross(w, ctr)
  return(w)
}

.minvolPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the minimum volatility portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)

  bounded <- !is.null(ctr$LB) || !is.null(ctr$UB)
  if (ctr$constraint[1] == "gross") {
    w <- .grossQP(Sigma, rep(0, n), ctr, n)
  } else if (!bounded) {
    tmp <- solve(Sigma, rep(1, n))
    w <- tmp/sum(tmp)
  } else {
    w <- .boundedQP(Sigma, rep(0, n), ctr, n)
  }
  w <- .finalizeWeights(w, ctr)
  .checkGross(w, ctr)
  return(w)
}

.invvolPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the inverse-volatility portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)

  ## The inverse-volatility portfolio is defined in closed form; it is always
  ## long-only and sums to one, so 'lo' and 'gross' hold automatically. Bounds
  ## do not enter its definition, but returning weights that violate the bounds
  ## the caller asked for is not an option either, so they are projected onto
  ## {LB <= w <= UB, w'1 = 1} and the caller is told.
  sig <- sqrt(diag(Sigma))
  if (any(sig <= 0)) {
    stop("the inverse-volatility portfolio requires strictly positive variances")
  }
  w <- 1/sig
  w <- w/sum(w)

  if (!is.null(ctr$LB) || !is.null(ctr$UB)) {
    wp <- .projectBox(w, ctr$LB, ctr$UB)
    if (max(abs(wp - w)) > 1e-12) {
      warning("the inverse-volatility weights violate 'LB'/'UB'; they have been ",
              "projected onto the bounds and no longer are exactly proportional ",
              "to the inverse volatilities")
    }
    w <- wp
  }
  return(w)
}

.ercExact <- function(Sigma, tol = 1e-10, maxit = 100000L) {
  ## Exact equal-risk-contribution weights. Solve y_i [Sigma y]_i = 1/N for
  ## y > 0 by cyclical coordinate descent -- each update is the positive root
  ## of Sigma_ii y_i^2 + (sum_{j != i} Sigma_ij y_j) y_i - 1/N = 0 -- and
  ## rescale to w = y / (y'1). This is the minimiser of the objective in
  ## '?optimalPortfolio' (it attains zero) and it is always strictly positive,
  ## so it satisfies the long-only constraint. See Spinu (2013).
  n <- dim(Sigma)[1]
  d <- diag(Sigma)
  if (any(!is.finite(d)) || any(d <= 0)) {
    stop("Sigma must have positive diagonal entries to compute the erc portfolio")
  }
  y <- 1/sqrt(d)
  ## Convergence is measured on the defining condition itself, relative to its
  ## 1/N target, and not by the size of the last step: the step is normalised
  ## by the largest component of y, so components orders of magnitude smaller
  ## were declared converged while still moving. Ill-conditioned matrices can
  ## need several thousand sweeps -- a fixed budget of a thousand returned
  ## risk contributions differing by tens of percentage points without saying so.
  resid <- Inf
  converged <- FALSE
  for (it in seq_len(maxit)) {
    for (i in 1:n) {
      b <- if (n > 1) sum(Sigma[i, -i] * y[-i]) else 0
      y[i] <- (-b + sqrt(b^2 + 4 * d[i]/n))/(2 * d[i])
    }
    resid <- max(abs(n * y * as.numeric(Sigma %*% y) - 1))
    if (is.finite(resid) && resid < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning("the equal-risk-contribution weights did not converge after ", maxit,
            " iterations; the risk contributions still depart from 1/N by up to ",
            format(resid, digits = 3), " in relative terms")
  }
  w <- as.numeric(y/sum(y))
  return(w)
}
.ercPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the equal-risk-contribution portfolio INPUTs
  ## Sigma : matrix (N x N) covariance matrix control : list of control
  ## parameters OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[2]
  ctr <- .ctrPortfolio(n, control)

  ## Without bounds, or under the long-only constraint, the exact solution is
  ## available in closed-ish form and is strictly positive. The squared
  ## deviation objective below is not convex: it has local minima on the
  ## boundary where a whole block of assets is dropped and the survivors are
  ## equalised at 1/(N - k), and slsqp started at the equally-weighted
  ## portfolio converges to one of them as soon as the correlation matrix has
  ## sizeable negative entries.
  if (ctr$constraint[1] == "none" || ctr$constraint[1] == "lo") {
    w <- .finalizeWeights(.ercExact(Sigma), ctr)
    .checkGross(w, ctr)
    return(w)
  }

  .pRC <- function(w) {
    Sigmaw <- crossprod(Sigma, w)
    pRC <- (w * Sigmaw)/as.numeric(crossprod(w, Sigmaw))
    d <- sum((pRC - 1/n)^2)
    return(d)
  }

  .gradERC <- function(w)
  {
    Sigmaw <- crossprod(Sigma, w)
    pRC <- (w * Sigmaw) / as.numeric(crossprod(w, Sigmaw))
    sig_p <- as.numeric(sqrt(crossprod(w, Sigmaw)))
    f <- pRC - 1/n
    g <- 2 * (sig_p^2 * (crossprod(Sigma, (w * f)) + f * Sigmaw) - 2 * Sigmaw * as.numeric(crossprod(w * f, Sigmaw))) / sig_p^4
  }

  ..grossContraint = NULL
  if (ctr$constraint[1] == "gross") {
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
  }

  ## the unconstrained erc portfolio is the natural extra starting value here
  starts <- .startingValues(ctr, Sigma)
  ## only a starting value here, so a non-converged run is not worth reporting
  starts <- unique(c(starts, list(.projectBox(suppressWarnings(.ercExact(Sigma)),
                                              ctr$LB, ctr$UB))))

  w <- .slsqpBest(starts = starts, fn = .pRC, gr = .gradERC,
                  hin = ..grossContraint, ctr = ctr)

  w <- .finalizeWeights(w, ctr)
  .checkGross(w, ctr)
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
  
  .gradMaxDiv <- function(w)
  {
    Sigmaw <- crossprod(Sigma, w)
    sig = sqrt(diag(Sigma))
    sig_p <- as.numeric(sqrt(crossprod(w, Sigmaw)))
    g <- (sig_p * sig - as.numeric(crossprod(w, sig)) * Sigmaw / sig_p) / sig_p^2
    g <- - g
  }
  
  ..grossContraint = NULL
  if (ctr$constraint[1] == "gross") {
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
  } else if (is.null(ctr$LB) && is.null(ctr$UB)) {
    ## nothing bounds the weights here, and the diversification ratio has no
    ## maximum when Sigma is singular
    .requirePD(Sigma, "the maximum-diversification ratio")
  }

  w <- .slsqpBest(starts = .startingValues(ctr, Sigma), fn = .divRatio,
                  gr = .gradMaxDiv, hin = ..grossContraint, ctr = ctr)

  w <- .finalizeWeights(w, ctr)
  .checkGross(w, ctr)
  return(w)
}

.riskeffPortfolio <- function(Sigma, semiDev, control = list()) {
  ## Compute the weight of the risk-efficient portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[2]
  ctr <- .ctrPortfolio(n, control)
  
  if (is.null(semiDev)) {
    stop("A vector of semideviations (semiDev) is required to compute the risk-efficient portfolio")
  }
  
  ## J is (N x 10): one column per decile of the semideviations, whatever the
  ## number of assets. Allocating it (N x N) instead made the loop below index
  ## past the last column, and hence fail, for every portfolio of fewer than
  ## ten assets. The lower edge of the first bin is -Inf so that every asset is
  ## assigned to exactly one decile, and deciles left empty by ties get a zero
  ## xi rather than an NA.
  pct <- as.numeric(quantile(semiDev, probs = seq(0.1, 1, 0.1)))
  lower <- c(-Inf, pct[1:9])
  epsilon <- vector("double", 10)
  J <- matrix(0, nrow = n, ncol = 10)

  for (i in 1:10) {
    pos <- semiDev > lower[i] & semiDev <= pct[i]
    J[pos, i] <- 1
    epsilon[i] <- if (any(pos)) median(semiDev[pos]) else 0
  }
  Jepsilon <- crossprod(t(J), epsilon)
  # DA Additional constraints used to stabilize optimization
  if (is.null(ctr$LB)) {
    ctr$LB <- (1/(2 * n)) * rep(1, n)
  }
  if (is.null(ctr$UB)) {
    ctr$UB <- (2/n) * rep(1, n)
  }
  ctr$w0 <- .projectBox(ctr$w0, ctr$LB, ctr$UB)

  .distRiskEff <- function(w) {
    Sigmaw <- crossprod(Sigma, w)
    d <- as.numeric(-crossprod(w, Jepsilon)/sqrt(crossprod(w, Sigmaw)))
    return(d)
  }
  
  .gradRiskEff <- function(w)
  {
    Sigmaw <- crossprod(Sigma, w)
    sig_p <- as.numeric(sqrt(crossprod(w, Sigmaw)))
    g <- (sig_p * Jepsilon - as.numeric(crossprod(w, Jepsilon)) * Sigmaw / sig_p) / sig_p^2
    g <- - g
  }
  
  ..grossContraint = NULL
  if (ctr$constraint[1] == "gross") {
    ..grossContraint = function(w) .grossConstraint(w, ctr$gross.c)
  }

  w <- .slsqpBest(starts = .startingValues(ctr, Sigma), fn = .distRiskEff,
                  gr = .gradRiskEff, hin = ..grossContraint, ctr = ctr)

  w <- .finalizeWeights(w, ctr)
  .checkGross(w, ctr)
  return(w)
}

.maxdecPortfolio <- function(Sigma, control = list()) {
  ## Compute the weight of the maximum decorrelation portfolio INPUTs Sigma :
  ## matrix (N x N) covariance matrix control : list of control parameters
  ## OUTPUTs w : vector (N x 1) weight
  n <- dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)
  Rho <- stats::cov2cor(Sigma)
  bounded <- !is.null(ctr$LB) || !is.null(ctr$UB)
  if (ctr$constraint[1] == "gross") {
    w <- .grossQP(Rho, rep(0, n), ctr, n)
  } else if (!bounded) {
    ## the maximum-decorrelation portfolio minimises w'Rw, not w'Sigma w:
    ## solving with Sigma here returned the minimum-variance portfolio
    tmp <- solve(Rho, rep(1, n))
    w <- tmp/sum(tmp)
  } else {
    w <- .boundedQP(Rho, rep(0, n), ctr, n)
  }
  w <- .finalizeWeights(w, ctr)
  .checkGross(w, ctr)
  return(w)
}

## Quadratic programs used by the convex portfolio problems.

.boundedQP <- function(Q, lin, ctr, n) {
  ## min 1/2 w'Q w - lin'w  s.t.  w'1 = 1 and whatever bounds were supplied.
  ## Every supplied bound enters the program: applying a bound afterwards, by
  ## projecting the solution of a differently-constrained problem onto the box,
  ## restores feasibility but does not minimise the objective.
  Amat <- matrix(rep(1, n), ncol = 1)
  bvec <- 1
  if (!is.null(ctr$LB)) {
    Amat <- cbind(Amat, diag(n))
    bvec <- c(bvec, ctr$LB)
  }
  if (!is.null(ctr$UB)) {
    Amat <- cbind(Amat, -diag(n))
    bvec <- c(bvec, -ctr$UB)
  }
  quadprog::solve.QP(Dmat = Q, dvec = lin, Amat = Amat, bvec = bvec, meq = 1)$solution
}

.grossQP <- function(Q, lin, ctr, n) {
  ## min 1/2 w'Q w - lin'w  s.t.  w'1 = 1, ||w||_1 <= gross.c and the bounds.
  ## Splitting w = p - n with p, n >= 0 turns the L1 budget into the linear
  ## constraint sum(p + n) <= gross.c; the relaxation is exact because any
  ## admissible w has a representation attaining it. The problem is therefore
  ## an ordinary quadratic program and is solved exactly, instead of being
  ## handed to a general-purpose nonlinear optimizer that has no attainable
  ## stopping criterion on it.
  Z <- cbind(diag(n), -diag(n))
  Dmat <- crossprod(Z, Q %*% Z)
  ## Z'QZ is singular by construction; a relative ridge makes it invertible for
  ## solve.QP and simultaneously favours the p_i n_i = 0 representation
  ridge <- 1e-09 * max(diag(Dmat))
  if (!is.finite(ridge) || ridge <= 0) {
    ridge <- 1e-09
  }
  Dmat <- Dmat + ridge * diag(2 * n)
  dvec <- as.numeric(crossprod(Z, lin))
  Amat <- cbind(c(rep(1, n), rep(-1, n)), -rep(1, 2 * n), diag(2 * n))
  bvec <- c(1, -ctr$gross.c, rep(0, 2 * n))
  if (!is.null(ctr$LB)) {
    Amat <- cbind(Amat, rbind(diag(n), -diag(n)))
    bvec <- c(bvec, ctr$LB)
  }
  if (!is.null(ctr$UB)) {
    Amat <- cbind(Amat, rbind(-diag(n), diag(n)))
    bvec <- c(bvec, -ctr$UB)
  }
  sol <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec,
                            meq = 1)$solution
  return(sol[1:n] - sol[(n + 1):(2 * n)])
}

## Constraints used by the optimizers. Both are stated in the convention
## required by nloptr >= 2.0.0, that is heq(w) == 0 and hin(w) <= 0.
.eqConstraint <- function(w) {
  return(sum(w) - 1)
}

.grossConstraint <- function(w, gross.c) {
  return(norm(as.matrix(w), type = "1") - gross.c)
}
