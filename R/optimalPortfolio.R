optimalPortfolio = function(Sigma, mu = NULL, semiDev = NULL, control = list()){ 
  #########################################################################################
  # Computes the weight of the chosen portfolio
  # INPUTs
  #   Sigma        : matrix (N x N) covariance matrix
  #   mu           : vector (N x 1) expected returns (optional)
  #   semiDev      : vector (N x 1) semi-deviation (optional)
  #   control      : a control list 
  #   The argument control is a list that can supply any of the following components
  #     type       : "mv", "minvol", "erc", "maxdiv", "riskeff"
  #     constraint : "none", "lo" (long only), "gross"
  #     gross.c    : default = 1.6 (130%-30%)
  #     gamma      : risk-averesion default = 0.89
  # OUTPUTs
  #   w            : vector (N x 1) optimal weights
  ######################################################################################### 
  if (missing(Sigma)){
    stop ('A covariance matrix (Sigma) is required')
  }
  if (!is.matrix(Sigma)){
    stop ('Sigma must be a matrix')
  }
  if (!isSymmetric(Sigma)){
    stop ('Sigma must be a symmetric matrix')
  }  

  ctr = .ctrPortfolio(control)
  
  if (ctr$type[1] == "mv"){
    w = .mvPortfolio(mu = mu, Sigma = Sigma, control = control)
  }
  else if (ctr$type[1] == "minvol"){
    w = .minvolPortfolio(Sigma = Sigma, control = control)
  }
  else if (ctr$type[1] == "erc"){
    w = .ercPortfolio(Sigma = Sigma, control = control)
  }
  else if (ctr$type[1] == "maxdiv"){
    w = .maxdivPortfolio(Sigma = Sigma, control = control)
  }
  else if (ctr$type[1] == "riskeff"){
    w = .riskeffPortfolio(Sigma = Sigma, semiDev = semiDev, control = control)
  }
  else if (ctr$type[1] == "invvol"){
    w = .invvolPortfolio(Sigma = Sigma, control = control)
  }
  else{
    stop('control$type is not well defined')
  }
  return(w)
}

.ctrPortfolio = function(control = list()){
  #########################################################################################  
  # Function used to control the list input
  # INPUTs
  #   control : a control list 
  #   The argument control is a list that can supply any of the following components
  #     type : "mv", "minvol", "erc", "maxdiv", "riskeff"
  #     constraint : "none", "lo", "gross"
  #     gross.c    : default = 1.6
  #     gamma      : default = 0.89
  # OUTPUTs
  #   control : list
  #########################################################################################  
  if (!is.list(control)){
    stop ('control must be a list') 
  }
  if (length(control) == 0){
    control = list(type = "mv", constraint = "none", gross.c = 1.6, 
                   LB = NULL, UB = NULL, w0 = NULL, gamma = 0.89)
  }
  nam = names(control)
  if (!("type" %in% nam) || is.null(control$type)){
    control$type = c("mv", "minvol", "erc", "maxdiv", "riskeff", "invvol")
  }
  if (!("constraint" %in% nam) || is.null(control$constraint)){
    control$constraint = c("none", "lo", "gross")
  }
  if (!("gross.c" %in% nam) || is.null(control$gross.c)){
    control$gross.c = 1.6
  }
  if (!("LB" %in% nam) || is.null(control$LB)){
    control$LB = NULL
  }
  if (!("UB" %in% nam) || is.null(control$UB)){
    control$UB = NULL
  }
  if (!("w0" %in% nam) || is.null(control$w0)){
    control$w0 = NULL
  }
  if (!("gamma" %in% nam) || is.null(control$gamma)){
    control$gamma = c(0.8773, 2.7063, 3.7950)
  }
  return(control)
}

.mvPortfolio = function(mu, Sigma, control = list()){
  #########################################################################################  
  # Compute the weight of the mean-variance portfolio
  # INPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #   control : list of control parameters
  # OUTPUTs
  #   w     : vector (N x 1) weight
  #########################################################################################    
  ctr = .ctrPortfolio(control)
  
  if (is.null(mu)){
    stop ('A vector of mean (mu) is required to compute the mean-variance portfolio')
  }
  if (ctr$constraint[1] == "none"){
    invSigmamu = solve(Sigma, mu)
    w = (1 / ctr$gamma[1]) * invSigmamu / sum(invSigmamu)  # David, rf = 0 !??
  }
  else if (ctr$constraint[1] == "lo"){
    n = dim(Sigma)[1]
    Dmat = ctr$gamma[1] * Sigma
    Amat = cbind(rep(1,n), diag(n))
    bvec = c(1, rep(0,n))
    w    = quadprog::solve.QP(Dmat = Dmat, dvec = mu, Amat = Amat, bvec = bvec, meq = 1)$solution   
  }
  else if (ctr$constraint[1] == "gross"){
    .meanvar = function(w){
      Sigmaw = crossprod(Sigma, w)
      opt    = -as.numeric( crossprod(mu, w) ) + 0.5 * ctr$gamma[1] * as.numeric( crossprod(w, Sigmaw) ) 
      return(opt)
    }
    
    n  = dim(Sigma)[1]
    w0 = ctr$w0
    if (is.null(w0)) {
      w0 = rep(1/n, n)  
    }
    w = nloptr::slsqp(x0 = w0, fn = .meanvar, hin = .grossConstraint, heq = .eqConstraint, nl.info = FALSE, 
              control = list(xtol_rel = 1e-18, check_derivatives = FALSE))$par
  }
  else{
    stop ('control$constraint not well defined')
  }
  return(w)
}

.minvolPortfolio = function(Sigma, control = list() ){
  #########################################################################################  
  # Compute the weight of the minimum volatility portfolio
  # INPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #   control : list of control parameters
  # OUTPUTs
  #   w     : vector (N x 1) weight
  #########################################################################################   
  ctr = .ctrPortfolio(control)
  n = dim(Sigma)[1]
  
  if (ctr$constraint[1] == "none"){
    tmp = solve(Sigma, rep(1, n))
    w   = tmp / sum(tmp)
  }
  else if (ctr$constraint[1] == "lo"){
    dvec = rep(0, n)
    Amat = cbind(rep(1, n), diag(n))
    bvec = c(1, rep(0, n))
    w    = quadprog::solve.QP(Dmat = Sigma, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)$solution    
  }
  else if (ctr$constraint[1] == "gross"){
    
    .minvol = function(w){
      Sigmaw = crossprod(Sigma, w)
      v = as.numeric(crossprod(w, Sigmaw))
      return(v)
    }
    
    n  = dim(Sigma)[1]
    w0 = ctr$w0
    if (is.null(w0)) {
      w0 = rep(1/n, n)  
    } 
    w = nloptr::slsqp(x0 = w0, fn = .minvol, hin = .grossConstraint, heq = .eqConstraint, 
              nl.info = FALSE, control = list(xtol_rel = 1e-18, check_derivatives = FALSE))$par
  }
  else{    
    stop ('control$constraint not well defined')    
  }
  return(w)  
}

.ercPortfolio = function(Sigma, control = list()){
  #########################################################################################  
  # Compute the weight of the equal-risk-contribution portfolio
  # INPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #   control : list of control parameters
  # OUTPUTs
  #   w     : vector (N x 1) weight
  #########################################################################################   
  ctr = .ctrPortfolio(control)
  
  n  = dim(Sigma)[2]
  w0 = ctr$w0
  if (is.null(w0)) {
    #w0 = rep(1/n, n)  
    w0 = 1 / sqrt(diag(Sigma))
    w0 = w0 / sum(w0)
  } 
  
  .pRC = function(w){
    Sigmaw = crossprod(Sigma, w)
    pRC    = (w * Sigmaw) / as.numeric( crossprod(w, Sigmaw) )
    d      = sum( (pRC - 1/n)^2 )
    return(d)
  }
  
  if (ctr$constraint[1] == "none"){
    w = nloptr::slsqp(x0 = w0, fn = .pRC, heq = .eqConstraint, lower = rep(0,n), 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE))$par  
  }
  else if(ctr$constraint[1] == "lo"){ 
    w = nloptr::slsqp(x0 = w0, fn = .pRC, heq = .eqConstraint, lower = rep(0,n), 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par  
  }
  else if(ctr$constraint[1] == "gross"){    
    w = nloptr::slsqp(x0 = w0, fn = .pRC, heq = .eqConstraint, lower = rep(0,n), 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par     
  }
  else{
    stop ('control$constraint not well defined')
  }
  return(w)
}

.maxdivPortfolio = function(Sigma, control = list()){ 
  #########################################################################################  
  # Compute the weight of weight of the maximum diversification portfolio
  # INPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #   control : list of control parameters
  # OUTPUTs
  #   w     : vector (N x 1) weight
  ######################################################################################### 
  ctr = .ctrPortfolio(control)
  
  n  = dim(Sigma)[2]
  w0 = ctr$w0
  if (is.null(w0)) {
    w0 = rep(1/n, n)  
  } 
  
  .divRatio = function(w){
    sig      = sqrt(diag(Sigma))
    Sigmaw   = crossprod(Sigma, w)
    divRatio = as.numeric(- crossprod(w, sig) / sqrt( crossprod(w, Sigmaw)) )
    return(divRatio)
  }
  if (ctr$constraint[1] == "none"){
    w = nloptr::slsqp(x0 = w0, fn = .divRatio, heq = .eqConstraint, 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par
  }
  else if(ctr$constraint[1] == "lo"){ 
    w = nloptr::slsqp(x0 = w0, fn = .divRatio, lower = rep(0,n), heq = .eqConstraint, 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par    
  }
  else if(ctr$constraint[1] == "gross"){
    w = nloptr::slsqp(x0 = w0, fn = .divRatio, hin = .grossConstraint, heq = .eqConstraint, 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par
  }
  else{
    stop ('control$constraint not well defined')
  }
  return(w)
}

.invvolPortfolio = function(Sigma, control = list()){
  #########################################################################################  
  # Compute the weight of the inverse-volatility portfolio
  # INPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #   control : list of control parameters
  # OUTPUTs
  #   w     : vector (N x 1) weight
  ######################################################################################### 
  
  sig = sqrt(diag(Sigma))
  w   = 1/sig
  w   = w / sum(w)
  return(w)
}

.riskeffPortfolio = function(Sigma, semiDev, control = list()){
  #########################################################################################  
  # Compute the weight of the risk-efficient portfolio
  # INPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #   control : list of control parameters
  # OUTPUTs
  #   w     : vector (N x 1) weight
  ######################################################################################### 
  ctr = .ctrPortfolio(control)
  
  if (is.null(semiDev)){
    stop ('A vector of semideviation (semiDev) is require to compute the risk-efficient portfolio')
  }
  
  n       = dim(Sigma)[2]
  pct     = c(0, quantile(semiDev, probs = seq(0.1, 1, 0.1))) 
  epsilon = vector("double", n)
  J       = matrix(rep(0, n^2), ncol = n)  
  
  for (i in 2:11) {
    pos = semiDev > pct[i-1] & semiDev <= pct[i]
    J[pos, i - 1] = 1
    epsilon[i - 1] = median(semiDev[pos])
  }
  Jepsilon = crossprod(t(J), epsilon)
  # !!! TOFIX !!! additional constraints used to stabilize optimization
  LB = (1 / (2 * n)) * rep(1, n)
  UB = (2 / n) * rep(1, n)
  
  w0 = ctr$w0
  if (is.null(w0)) {
    w0 = (UB - LB)
    w0 = w0 / sum(w0)  
  } 
  
  .distRiskEff = function(w){
    Sigmaw = crossprod(Sigma, w)
    d      = as.numeric(-crossprod(w, Jepsilon) / sqrt( crossprod(w, Sigmaw)))
    return(d)
  }
  
  .eqConstraintReff = function(w){
    return(sum(w))
  }
  
  if (ctr$constraint[1] == "none"){
    w = nloptr::slsqp(x0 = w0, fn = .distRiskEff, heq = .eqConstraint, lower = LB, upper = UB, 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par
  }
  else if(ctr$constraint[1] == "lo"){ 
    w = nloptr::slsqp(x0 = w0, fn = .distRiskEff, heq = .eqConstraint, lower = LB, upper = UB, 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par    
  }
  else if(ctr$constraint[1] == "gross"){
    w = nloptr::slsqp(x0 = w0, fn = .distRiskEff, heq = .eqConstraint, lower = LB, upper = UB, 
                      nl.info = FALSE, control = list(xtol_rel = 1e-8, check_derivatives = FALSE, maxeval = 2000))$par
  }
  else{
    stop ('control$constraint not well defined')
  }
  return(w)
}

#########################################################################################  
## Constraints used by the optimizers
#########################################################################################  
.eqConstraint = function(w){
  return(sum(w) - 1)
}

# DA here 1.6 is hard-coded, this should be changed
.grossConstraint = function(w){
  return(1.6 - norm(as.matrix(w), type = "1") )
}

