meanEstimation = function(rets, control = list()){
  #########################################################################################  
  # Compute the estimation of the mean 
  # INPUTs
  #   rets     : matrix (T x N) returns
  #   control  : a control list 
  #   The argument control is a list that can supply any of the following components
  #     type   : "naive", "ewma", "bs", "mart"  default = "naive"
  #     lambda : default = 0.94 
  # OUTPUTs
  #   mu : vector (N x 1) mean
  #########################################################################################
  
  if (missing(rets)){
    stop ('rets is missing')
  }
  if (!is.matrix(rets)){
    stop ('rets must be a (T x N) matrix')
  }
    
  ctr = .ctrMean(control)
  
  if (ctr$type[1] == "naive"){
    mu = .naiveMean(rets = rets)
  }
  else if (ctr$type[1] == "ewma"){
    mu = .ewmaMean(rets = rets, lambda = ctr$lambda)
  }
  else if (ctr$type[1] == "bs"){
    mu = .bsMean(rets = rets)
  }
  else if (ctr$type[1] == "mart"){
    mu = .martMean(rets = rets)
  }
  else{
    stop('control$type is not well defined')
  }
  return(mu)
} 

.ctrMean = function(control = list()){
  #########################################################################################  
  # Function used to control the list input
  # INPUTS
  #   control : a control list  
  #   The argument control is a list that can supply any of the following components
  #     type   : default = "naive"
  #     lambda : default = 0.94
  # OUTPUTs
  #   control : a list
  #########################################################################################  
  if (!is.list(control)){
    stop ('control must be a list') 
  }
  if (length(control) == 0){
    control = list(type = "naive", lambda = 0.94)
  }
  nam = names(control)
  if (!("type" %in% nam) || is.null(control$type)){
    control$type = c("naive", "ewma", "bs", "mart")
  }
  if (!("lambda" %in% nam) || is.null(control$lambda)){
    control$lambda = 0.94
  }
  return(control) 
}

.naiveMean = function(rets){
  #########################################################################################  
  # Compute the naive mean 
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   mu : vector (N x 1) mean
  #########################################################################################  
  mu = colMeans(rets)
  return(mu)
}

.ewmaMean = function(rets, lambda){
  #########################################################################################  
  # Compute the exponential weighted moving average 
  # INPUTs
  #   rets   : matrix (T x N) returns
  #   lambda : scalar
  # OUTPUTs
  #   mu : vector (N x 1) mean
  # NOTE
  #   the dates should be in ascending order, i.e. ordest at the top and newest at the bottom
  #########################################################################################  
  t  = dim(rets)[1]
  w  = lambda^(t:1)
  w  = w / sum(w)
  w  = matrix(data = w, nrow = t, ncol = dim(rets)[2], byrow = FALSE)
  mu = colSums(w * rets)
  return(mu)
}

.bsMean = function(rets){
  #########################################################################################  
  # Compute the Bayes-Stein mean
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   mu : vector (N x 1) mean
  #########################################################################################  
  # !!! TOFIX !!! estime pour l'instant une matrice de var-cov de facon standard
  # verifier avec david si on laisse le choix en input
  mu    = colMeans(rets)
  Sigma = cov(rets)
  invSigma = solve(Sigma)
  N = dim(Sigma)[1]
  T = length(mu)
  i = rep(1,N)
  invSigmai = crossprod(invSigma, i)
  w_min = (invSigmai) / as.numeric(crossprod(i, invSigmai))
  mu_min = crossprod(mu, w_min)
  invSigmaMu = crossprod(invSigma, mu - mu_min)
  phi = (N+2) / ((N+2) + T * crossprod(mu - mu_min, invSigmaMu ))
  phi = max(min(phi, 1), 0)
  mu = (1 - phi) * mu + phi * mu_min
  return(mu)
}

.martMean = function(rets){
  #########################################################################################  
  # Compute the Martinelli mean
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   mu   : vector (N x 1) mean
  #########################################################################################  
  mu = apply(rets, 2, sd)
  return(mu)
}