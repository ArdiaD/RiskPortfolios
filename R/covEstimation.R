covEstimation = function(rets, control = list()){
  #########################################################################################  
  # Compute the chosen covariance matrix
  # INPUTs
  #   rets : matrix (T x N) returns
  #   control : a control list 
  #   The argument control is a list that can supply any of the following components
  #     type   : "naive", "ewma", "lw", "factor", "rtm", "const", "cor", "oneparm", "diag", "large"
  #     lambda : 
  #     K      : 
  # OUTPUTs
  #   Sigma : matrix (N x N) variance-covariance 
  #########################################################################################   
  
  if (missing(rets)){
    stop ('rets is missing')
  }
  if (!is.matrix(rets)){
    stop ('rets must be a (T x N) matrix')
  }
    
  ctr = .ctrCov(control)
  
  if (ctr$type[1] == "naive"){
    Sigma = .naiveCov(rets = rets)
  }
  else if (ctr$type[1] == "ewma"){
    Sigma = .ewmaCov(rets = rets, lambda = ctr$lambda)
  }
  else if (ctr$type[1] == "lw"){
    Sigma = .lwCov(rets = rets)
  }
  else if (ctr$type[1] == "factor"){
    Sigma = .factorCov(rets = rets, K = ctr$K)
  }
  else if (ctr$type[1] == "const"){
    Sigma = .constCov(rets = rets)
  }
  else if (ctr$type[1] == "cor"){
    Sigma = .corCov(rets = rets)
  }  
  else if (ctr$type[1] == "oneparm"){
    Sigma = .oneparmCov(rets = rets)
  }  
  else if (ctr$type[1] == "diag"){
    Sigma = .diagCov(rets = rets)
  }
  else if (ctr$type[1] == "large"){
    Sigma = .largeCov(rets = rets)
  }
  else if (ctr$type[1] == "bs"){
    Sigma = .bsCov(rets = rets)
  }
  else{    
    stop('control$type is not well defined')
  }
  return(Sigma)
}  

.ctrCov = function(control = list()){
  #########################################################################################  
  # Function used to control the list input
  # INPUTS
  #   control  : a control list  
  #   The argument control is a list that can supply any of the following components
  #     type : default = "naive"
  #     lambda : default = 0.94
  #     K      : default = 1
  # OUTPUTs
  #   control : a list
  #########################################################################################
  if (!is.list(control)){
    stop ('control must be a list') 
  }
  if (length(control) == 0){
    control = list(type = "naive", lambda = 0.94, K = 1)
  }
  nam = names(control)
  if (!("type" %in% nam) || is.null(control$type)){
    control$type = c("naive", "ewma", "lw", "factor", "rtm", "const", "cor", "oneparm", "diag", "large", "bs")
  }
  if (!("lambda" %in% nam) || is.null(control$lambda)){
    control$lambda = 0.94
  }
  if (!("K" %in% nam) || is.null(control$K)){
    control$K = 1
  }
  return(control)  
}

.naiveCov = function(rets){
  #########################################################################################  
  # Compute the naive covariance matrix
  # INPUTs
  #   rets  : matrix (T x N) returns
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #########################################################################################    
  # DA here we could check that the dimension doesn't lead to bullshit outputs
  Sigma = cov(rets)
  return(Sigma)
}

.ewmaCov = function(rets, lambda){
  #########################################################################################  
  # Compute the exponential weighted moving average covariance matrix
  # INPUTs
  #   rets : matrix (T x N) returns
  #   lambda : 
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #########################################################################################  
  # DA check speed here (use just-in-time compilation or compiler package)
  
#   # Infinite sample
#   t = nrow(rets)
#   n = ncol(rets)
#   Sigma = matrix(0, ncol = n, nrow = n)
##   Sigma = cov(rets)

#   for (i in 1 : t){
#     Sigma  = lambda * Sigma + (1 - lambda) * outer( as.double(rets[i, ]), as.double(rets[i, ]) )
#   }
#   
  ## Finite sample ewma
  t     = nrow(rets)
  Sigma = cov(rets)
#   Sigma = matrix(0, ncol = dim(rets)[2], dim(rets)[2])
  mu    = colMeans(rets)
  shiftRets = sweep(rets, 2, mu, '-')
  for (i in 1 : t){
    r = as.double(shiftRets[i, ])
    r2 = outer(r,  r) 
    Sigma  = (1 - lambda) / (1-lambda ^ t) * r2  + lambda * Sigma
  }
  return(Sigma)  
}

.constCov = function(rets){
  ######################################################################################  
  #  Shrinkage of the covariance matrix using constant correlation matrix
  # INPUTs
  #   rets : matrix (T x N) returns 
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #######################################################################################  
  n      = dim(rets)[2]
  tmpMat = matrix(rep(1, n^2), ncol = n)
  
  rho     = mean(cor(rets)[lower.tri(tmpMat, diag = FALSE)] )
  R       = rho * tmpMat
  diag(R) = 1
  
  std     = apply(rets, 2, sd)
  diagStd = diag(std)
  Sigma   = diagStd %*% R %*% diagStd
  
  return(Sigma) 
}

.factorCov = function(rets, K){
  #########################################################################################  
  #  Compute the covariance matrix using K-factor approach
  # INPUTs
  #   rets : matrix (T x N) returns
  #   K : [scalar] number of factors
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  # NOTE
  #   Matlab function is factoran (statistical toolbox)
  #########################################################################################  
  # !!! TOFIX !!! valider avec Matlab
  std   = apply(rets, 2, sd)
  sigma = cov(rets)
  
  loading      = factanal(rets, K)$loadings
  uniquenesses = factanal(rets,K)$uniquenesses
  
  R       = tcrossprod(loading) + diag(uniquenesses) 
  diagStd = diag(std)
  Sigma   = diagStd %*% R %*% diagStd
  
  return(Sigma)
}

.lwCovElement = function(rets, type){
  #########################################################################################  
  #  Computes the common elements of the functions for Ledoit-Wolf
  # INPUTs
  #   rets : matrix (T x N) returns
  #   type : 'large', 'lw' or else
  # OUTPUTs
  #   list : useful outputs
  #########################################################################################  
  t         = dim(rets)[1]
  n         = dim(rets)[2]
  mu        = colMeans(rets)
  shiftRets = sweep(rets, 2, mu, "-")
  y         = shiftRets^2
  
  if (type == "large" ||type == "lw"){
#     mkt     = rowMeans(rets)
    mkt     = rowMeans(shiftRets)
      
    smple   = cov( cbind(rets, mkt) ) * (t - 1) / t
    covmkt  = smple[1 : n, n + 1]
    covmkt_ = matrix( rep(covmkt, n) , ncol = n, byrow = FALSE)
    varmkt  = as.numeric(smple[n + 1, n + 1])
    smple   = smple[- (n + 1), - (n + 1)]
    
    prior       = outer(covmkt, covmkt) / varmkt
    diag(prior) = diag(smple)
    
    lwCovElement = list(t = t, n = n, mu = mu, shiftRets = shiftRets, y = y, 
                        smple = smple, mkt = mkt, covmkt = covmkt, 
                        covmkt_ = covmkt_, varmkt = varmkt, prior = prior)
  }
  else{
    smple = (1 / t) * crossprod(shiftRets)
    
    lwCovElement = list(t = t, n = n, mu = mu, shiftRets = shiftRets, y = y, 
                        smple = smple, mkt = NULL, covmkt = NULL, 
                        covmkt_ = NULL, varmkt = NULL, prior = NULL)
  }
  return(lwCovElement)
}

.lwCov = function(rets){
  #########################################################################################  
  #  Shrinkage of the covariance matrix towards market
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #########################################################################################  
  lwCovElement = .lwCovElement(rets, type = "lw")
  
  t         = lwCovElement$t
  n         = lwCovElement$n
  mu        = lwCovElement$mu
  shiftRets = lwCovElement$shiftRets
  mkt       = lwCovElement$mkt
  covmkt    = lwCovElement$covmkt
  varmkt    = lwCovElement$varmkt
  smple     = lwCovElement$smple
  prior     = lwCovElement$prior
  y         = lwCovElement$y
  z         = sweep(shiftRets, 1, mkt, "*")
   
  # Phi hat
  phiMat  = crossprod(y) / t - 2 * crossprod(shiftRets) * smple / t + smple^2
  phi     = sum(apply(phiMat, 2, sum))
 
  # Rho hat  
  rhoMat1       = 1 / t * sweep(crossprod(y, z), 2, covmkt, "*") / varmkt  
  rhoMat3       = 1 / t * crossprod(z) * outer(covmkt, covmkt) / varmkt^2
  rhoMat        = 2 * rhoMat1 - rhoMat3 - prior*smple
  diag(rhoMat)  = diag(phiMat)
  rho           = sum(apply(rhoMat, 2, sum))

  # Gamma hat 
  gamma = norm(smple - prior, "F") ^ 2
  
  # Kappa hat
  kappa = (phi - rho) / gamma
  
  # shrinkage value
  shrinkage = pmax(0, pmin(1, kappa / t))
  
  # Sigma hat
  Sigma = shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}
  
.largeCov = function(rets){
  #########################################################################################  
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  # !!! TOFIX !!! largeCov ou puis-je trouver les sources pour formules ?
  #########################################################################################    
  lwCovElement = .lwCovElement(rets, type = "large")
  
  t         = lwCovElement$t
  n         = lwCovElement$n
  mu        = lwCovElement$mu
  shiftRets = lwCovElement$shiftRets
  
  y = lwCovElement$y
  
  mkt     = lwCovElement$mkt
  covmkt  = lwCovElement$covmkt
  covmkt_ = lwCovElement$covmkt_
  varmkt  = lwCovElement$varmkt
  smple   = lwCovElement$smple
  prior   = lwCovElement$prior
  
  z = sweep(shiftRets, 1, mkt, "*")
  
  d = 1 / n * norm(smple - prior, "F")^2
  r2 = 1 / n / t^2 * sum( apply(crossprod(y, y), 2, sum)) - 1 / n / t * sum(apply(smple^2, 2, sum))
  
  phidiag = 1 / n / t^2 * sum(apply(y^2, 2, sum)) - 1 / n / t * sum(diag(smple)^2)
 
  v1 = 1 / t^2 * crossprod(y,z) - 1 / t *  covmkt_ * smple
  
  phioff1 = 1 / n * sum(apply(v1 * t(covmkt_), 2, sum)) / varmkt - 1 / n * sum(diag(v1) *covmkt) / varmkt
  
  v3 = 1 / t^2 * crossprod(z,z) - 1 / t * varmkt * smple
  
  phioff3 = 1 / n * sum(apply(v3 * tcrossprod(covmkt, covmkt), 2, sum)) / varmkt^2 - 1 / n * sum(diag(v3) * covmkt^2) / varmkt^2
  
  phioff = 2 * phioff1 - phioff3
  phi = phidiag + phioff
  
  Sigma = (r2 - phi) / d * prior + (1 - (r2 - phi) / d) * smple
  
  return(Sigma)
}

.corCov = function(rets){
  #########################################################################################  
  #  Shrinkage of the covariance matrix towards constant correlation matrix by Ledoit-Wolf
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #########################################################################################  
  lwCovElement = .lwCovElement(rets, type = "cor")
  
  t             = lwCovElement$t
  n             = lwCovElement$n
  mu            = lwCovElement$mu
  shiftRets     = lwCovElement$shiftRets
  y             = lwCovElement$y
  smple         = lwCovElement$smple
  var           = diag(smple)
  sqrtvar       = sqrt(var)
  outerSqrtVar  = outer(sqrtvar, sqrtvar)
  
  rBar          = (sum(sum(smple / outerSqrtVar)) - n) / (n * (n - 1)) 
  prior         = rBar * outerSqrtVar
  diag(prior)   = var
  
  # phi hat
  phiMat  =crossprod(y) / t - 2 * crossprod(shiftRets) * smple / t + smple^2
  phi     = sum(apply(phiMat, 2, sum))
  
  # rho hat

  term1         = crossprod(shiftRets^3,shiftRets) / t
  term2         = sweep(smple, 1, diag(smple), "*")
  #term3         = sweep(smple, 1, var, "*") # lorsqu'on developpe on pourrait annuler les term3 et term4 ... R
  #term4         = sweep(smple, 1, var, "*")
  rhoMat        = term1 - term2 #- term3 + term4
  diag(rhoMat)  = 0
  rho           = sum(diag(phiMat)) + rBar * sum(sum( outer(1 / sqrtvar, sqrtvar) * rhoMat)) 

  # gamma hat
  gamma = norm(smple - prior, type = "F")^2
  
  # kappa hat
  kappa = (phi - rho) / gamma
  
  # shrinkage value
  shrinkage = pmax(0, pmin(1, kappa / t))
  
  # Sigma hat
  Sigma = shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}


.diagCov = function(rets){
  #########################################################################################  
  #  Shrinks towards diagonal matrix
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #########################################################################################  
  lwCovElement = .lwCovElement(rets, type = "diag")
  
  t         = lwCovElement$t
  n         = lwCovElement$n
  mu        = lwCovElement$mu
  shiftRets = lwCovElement$shiftRets
  y         = lwCovElement$y
  smple     = lwCovElement$smple
  prior     = diag(diag(smple))
  
  # phi hat
  phiMat  = crossprod(y) / t - 2 * (crossprod(shiftRets)) * smple / t + smple^2
  phi     = sum(apply(phiMat, 2, sum))
  
  # rho hat
  rho = sum(diag(phiMat))

  # gamma hat
  gamma = norm(smple - prior, "F")^2
  
  # kappa hat
  kappa = (phi - rho) / gamma
  
  # shrinkage value
  shrinkage = pmax(0, pmin(1, kappa / t))
  
  # Sigma hat
  Sigma = shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}

.oneparmCov = function(rets){
  #########################################################################################  
  #  Shrinks towards one-parameter matrix
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  #   Sigma : matrix (N x N) covariance matrix
  #########################################################################################  
  lwCovElement = .lwCovElement(rets, type = "oneparm")
  
  t         = lwCovElement$t
  n         = lwCovElement$n
  mu        = lwCovElement$mu
  shiftRets = lwCovElement$shiftRets
  smple     = lwCovElement$smple
  y         = lwCovElement$y
  meanvar   = mean(diag(smple))
  prior     = meanvar * diag(n)

  # phi hat
  phiMat  = crossprod(y) / t - 2 * (crossprod(shiftRets)) * smple / t + smple^2
  phi     = sum(apply(phiMat, 2, sum))
  
  # gamma hat
  gamma = norm(smple - prior,  type = "F")^2
  
  # kappa hat
  kappa = phi / gamma
  
  # shrinkage value
  shrinkage = pmax(0, pmin(1, kappa / t))
  
  # Sigma hat
  Sigma = shrinkage * prior + (1 - shrinkage) * smple
  
  return(Sigma)
}


.bsCov = function(rets){
  #########################################################################################  
  # Compute the Bayes-Stein covariance matrix
  # INPUTs
  #   rets : matrix (T x N) returns
  # OUTPUTs
  # Sigma : Matrix (N x N) covariance
  # David : see Kolusheva, Daniela. (July 2008) Out-of-sample Performance of Asset Allocation Strategies
  #########################################################################################  
  lwCovElement = .lwCovElement(rets, type = "lw")
  t         = lwCovElement$t
  n         = lwCovElement$n

  mu    = colMeans(rets)
  Sigma = cov(rets)
  invSigma = solve(Sigma)

  i = rep(1,n)
  invSigmai = crossprod(invSigma, i)
  w_min = (invSigmai) / as.numeric(crossprod(i, invSigmai))
  mu_min = crossprod(mu, w_min)
  invSigmaMu = crossprod(invSigma, mu - mu_min)
  phi = (n+2) / ((n+2) + t * crossprod(mu - mu_min, invSigmaMu ))
  phi = max(min(phi, 1), 0)
  
  tau = t * phi / (1 - phi)
  
  Sigma = Sigma*(1 + 1 / (t + tau)) + tau / (t * (t + 1 + tau)) * outer(i,i) /as.numeric( crossprod(i, invSigmai))
  
  return(Sigma)
}