semidevEstimation = function(rets, control = list()){
  #########################################################################################  
  # Compute the chosen semi-deviation
  # INPUTs
  #   rets     : matrix (T x N) returns
  #   control  : a control list 
  #   The argument control is a list that can supply any of the following components
  #     type   : "naive", "ewma"
  #     lambda : default = 0.94 
  # OUTPUTs
  #   semiDev  : vector (N x 1) semi-deviations
  #########################################################################################   
  
  if (missing(rets)){
    stop ('rets is missing')
  }
  if (!is.matrix(rets)){
    stop ('rets must be a (T x N) matrix')
  }
  
  ctr = .ctrSemidev(control)
  
  if (ctr$type[1] == "naive"){
    semiDev = .naiveSemiDev(rets = rets)
  }
  else if (ctr$type[1] == "ewma"){
    semiDev = .ewmaSemiDev(rets = rets, lambda = ctr$lambda)
  }
  else{
    stop('control$type is not well defined')
  }
  return(semiDev)
}

.ctrSemidev = function(control = list()){
  
  if (!is.list(control)){
    stop ('control must be a list') 
  }
  if (length(control) == 0){
    control = list(type = "naive", lambda = 0.94)
  }
  nam = names(control)
  if (!("type" %in% nam) || is.null(control$type)){
    control$type = c("naive", "ewma")
  }
  if (!("lambda" %in% nam) || is.null(control$lambda)){
    control$lambda = 0.94
  }
  return(control)
}

.naiveSemiDev = function(rets){
  #########################################################################################  
  # Compute the naive semi-deviation
  # INPUTs
  #   rets    : matrix (T x N) returns
  # OUTPUTs
  #   semiDev : vector (N x 1) semi-deviation
  #########################################################################################    
  semiDev = .ewmaSemiDev(rets, lambda = 1)
  return(semiDev)
}

.ewmaSemiDev = function(rets, lambda){
  #########################################################################################  
  # Compute the ewma semi-deviation
  # INPUTs
  #   rets    : matrix (T x N) returns
  # OUTPUTs
  #   semiDev : vector (N x 1) semi-deviation
  #########################################################################################    
  t       = dim(rets)[1]
  n       = dim(rets)[2]
  semiDev = vector("double", n)
  mu      = colMeans(rets)
  w       = lambda^(t:1) 
  for (j in 1 : n){
    retsj  = rets[,j]
    muj    = mu[j]
    idx    = retsj < muj
    wj     = w[idx] / sum(w[idx])
    semiDev[j] = sqrt(sum(wj * (retsj[idx] - muj)^2))
  }
  return(semiDev)
}