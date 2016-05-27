impliedReturns = function(rets = NULL, mu = NULL, Sigma = NULL, semiDev = NULL, w = NULL, control = list()){
  #########################################################################################  
  # Compute the implied expected returns 
  #########################################################################################
  
  ctr = .ctrImpliedReturns(control)
  
  if (is.null(mu)){
    mu = meanEstimation(rets = rets, control = ctr$control.mean)
  }
  if (is.null(Sigma) & is.null(semiDev)){
    Sigma = covEstimation(rets = rets, control = ctr$control.cov)
  }
  if (is.null(semiDev) & is.null(Sigma)){
    semiDev = semidevEstimation(rets, control = ctr$control.semidev)
  }
  if (is.null(w)){
    w = optimalPortfolio(mu = mu, Sigma = Sigma, semiDev = semiDev, control = ctr$control.ptf)
  }
 
  if (ctr$type[1] == "bl"){
    imu = .blImpliedReturns(Sigma = Sigma, w = w, gamma = ctr$gamma[1])
  } 
  else if (ctr$type[1] == "regression"){
    if (ctr$reg[1] == "standard"){
      imu = .regImpliedReturns(mu = mu, Sigma = Sigma, w = w)
    }
    else if (ctr$reg[1] == "robust"){
      imu = .robregImpliedReturns(mu = mu, Sigma = Sigma, w = w)
    }
    else if (ctr$reg[1] == "constraint"){
      imu = .constraintImpliedReturns(mu = mu, Sigma = Sigma, w = w)
    }
    else{
      stop('control$reg not well defined')
    }
  }
  else{
    stop('control$type not well defined')
  }
  
  out = list(imu = imu, mu = mu, Sigma = Sigma, w = w)
  return(out) 
}

.ctrImpliedReturns = function(control = list()){
  
  if (!is.list(control)){
    stop ('control must be a list') 
  }
  if (length(control) == 0){
    control = list(type = "regression", reg = "standard", control.mean = NULL, control.cov = NULL, control.ptf = NULL)
  }
  nam = names(control)
  if (!("type" %in% nam) || is.null(control$type)){
    control$type = c("regression", "bl")
  }
  if (!("gamma" %in% nam) || is.null(control$gamma)){
    control$gamma = .ctrPortfolio()$gamma
  }
  if (!("reg" %in% nam) || is.null(control$reg)){
    control$reg = c("standard", "robust", "constraint")
  }
  if (!("control.mean" %in% nam) || is.null(control$control.mean)){
    control$control.mean = .ctrMean()
  }
  if (!("control.cov" %in% nam) || is.null(control$control.cov)){
    control$control.cov = .ctrCov()
  }
  if (!("control.semidev" %in% nam) || is.null(control$control.semidev)){
    control$control.semidev = .ctrSemidev()
  }
  if (!("control.ptf" %in% nam) || is.null(control$control.ptf)){
    control$control.ptf = .ctrPortfolio()
  }
  return(control)
}

.blImpliedReturns = function(Sigma, w, gamma){
  
  imu = gamma * crossprod(Sigma, w)
  return(imu)
}

.regImpliedReturns = function(mu, Sigma, w){
  
  # DA this is unconstrained for the moment
  imu = as.vector(lm(mu ~ (Sigma * w))$fitted) 
  return(imu)
}

.robregImpliedReturns = function(mu, Sigma, w){
  
  # DA robust regression here
  imu = as.vector(MASS::rlm(mu ~ (Sigma * w))$fitted) 
  return(imu)
}


.constraintImpliedReturns = function(mu, Sigma, w){

  x      = crossprod(Sigma, w)
  xBar   = mean(x)
  shiftX = x - xBar
  yBar   = mean(mu)
  shiftY = mu - yBar
  b      = max( crossprod(shiftX, shiftY) / sum(x^2), 0)
  a      = yBar - b * xBar
  muT    = a + b * x
  return(muT)  
}
