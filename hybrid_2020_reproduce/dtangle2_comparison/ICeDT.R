
meanFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = mean)
  return(out)
}

varFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = var)
  return(out)
}

#-------------------------------------------------------------------#
# initial estimate of cell type composition
# yr = c(tumor purity, expression of nG genes)
# Z is the expression signature matrix
#-------------------------------------------------------------------#

lmInit<-function(yr, Z){
  rho_t = yr[1]
  y     = yr[-c(1)]
  
  initMod = lm(y~Z)
  
  rho  = coef(initMod)[-1]
  wneg = which(rho < 0.005)
  rho[wneg] = 0.05
  rho = (rho/sum(rho))*(1-rho_t)
  
  return(rho)
}

#-------------------------------------------------------------------#
#  initial estimate of variance for consistent or abberant genes
#-------------------------------------------------------------------#

sigmaInit <- function(x, Z, nG, var_wgt){
  nCT    = ncol(Z)
  Y      = x[c(1:nG)]
  logY   = log(Y)
  rho    = x[c((nG+1):(nG+nCT))]
  eta_ij = Z %*% rho
  
  resid = Y - eta_ij
  Q3val = quantile(abs(resid), probs = c(0.75),na.rm=TRUE)
  g2use = which(abs(resid) > Q3val)
  
  init_val = c(0,0)
  
  # g2use is the initial set of abberant genes. It represents the first 
  # look at aberrant genes by selecting them from the the upper 
  # quartile of residuals from initial fit (worst fit 25%). 
  # The remaining 75% are considered "consistent marker genes".
  
  log_eta_ij = log(eta_ij)
  init_val[1] = sigma2_Update(logY = logY[-c(g2use)], 
                              log_eta_ij = log_eta_ij[-c(g2use)],
                              EM_wgt = rep(1,(nG-length(g2use))), 
                              AB_Up = FALSE, var_wgt=var_wgt)
  
  init_val[2] = sigma2_Update(logY = logY[c(g2use)], 
                              log_eta_ij = log_eta_ij[c(g2use)], 
                              EM_wgt = rep(0,length(g2use)), 
                              AB_Up = TRUE, var_wgt=var_wgt)
  return(init_val)
}


#-------------------------- Sigma Updates --------------------------#

sigma2_obj_function <- function(x, logY, log_eta_ij, var_wgt, EM_wgt, 
                                AB_Up=FALSE){
  sigma2_wgt  = var_wgt*x
  
  mu_ij = log_eta_ij - sigma2_wgt/2
  
  if(AB_Up==FALSE){
    out = sum(EM_wgt*dnorm(logY, mean=mu_ij, sd=sqrt(sigma2_wgt),log=TRUE))
  } else {
    out = sum((1-EM_wgt)*dnorm(logY,mean=mu_ij, sd=sqrt(sigma2_wgt), log=TRUE))
  }
  return(out)
}

sigma2_Update <- function(logY, log_eta_ij, EM_wgt, AB_Up=FALSE, 
                          var_wgt=NULL){
  
  if(is.null(var_wgt)){
    varTheta_ij = logY - log_eta_ij
    
    if(AB_Up==FALSE){
      Sigma2_Up = 2*(sqrt((sum(EM_wgt*(varTheta_ij^2))/sum(EM_wgt))+1)-1)
    } else {
      Sigma2_Up = 2*(sqrt((sum((1-EM_wgt)*(varTheta_ij^2))/sum(1-EM_wgt))+1)-1)
    }
  }else{
    optimizeOut = optimize(f = sigma2_obj_function, logY = logY,
                           log_eta_ij=log_eta_ij, var_wgt=var_wgt, 
                           AB_Up=AB_Up, EM_wgt=EM_wgt, 
                           lower=0, upper=10, maximum = TRUE)
    
    Sigma2_Up = optimizeOut$maximum
  }
  
  return(Sigma2_Up)
}

#------------------------ Support Functions ------------------------#

correctRho <- function(est, total){
  rho  = est
  wneg = which(rho<0.005)
  rho[wneg] = 0.005
  
  rho = total*rho/sum(rho)
  
  return(rho)
}

correctRho_v2 <- function(est){
  rho = est
  wneg = which(rho<0.005)
  rho[wneg] = 0.005
  
  return(rho)
}

extractFunction <- function(compList, element){
  out = mapply(compList, FUN = function(x){ get(element, x) })
  return(out)
}

#----------------------- Gradient Functions ------------------------#

gradFunc_noPurity <- function(x, logY, Z, sigma2C, sigma2A, EM_wgt){
  eta_ij = Z %*% x
  mu_ijC = log(eta_ij) - sigma2C/2
  d_ijC  = logY - mu_ijC
  
  mu_ijA = log(eta_ij) - sigma2A/2
  d_ijA  = logY - mu_ijA
  
  c1 = c(d_ijC*EM_wgt/(sigma2C*eta_ij)) + c(d_ijA*(1-EM_wgt)/(sigma2A*eta_ij))
  
  out = t(Z) %*% matrix(c1,ncol=1)
  
  return(out)
}

gradFunc_givenPurity <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho = c(x, 1 - rho_i0 -sum(x))
  eta_ij = Z%*%rho
  
  mu_ijC = log(eta_ij) - sigma2C/2
  d_ijC  = logY - mu_ijC
  
  mu_ijA = log(eta_ij) - sigma2A/2
  d_ijA  = logY - mu_ijA
  
  c1 = c(d_ijC*EM_wgt/(sigma2C*eta_ij)) + c(d_ijA*(1-EM_wgt)/(sigma2A*eta_ij))
  
  # the derivative was derived using all cell type composition of non-tumor 
  # cell types. Since the proportion of last cell type is 
  # rho_K = constant - sum(x), its derivative is the original derivative
  # multiples d [rho_K] / d [rho_k], for k = 1 to K-1, which is - 1 
  # so it is equivalent to replace Z by Z_star
  nCT    = ncol(Z)
  Z_star = Z %*% rbind(diag(rep(1,(nCT-1))), rep(-1,(nCT-1)))
  
  out = t(Z_star) %*% matrix(c1,ncol=1)
  
  return(out)
}

#----------------- Likelihood (no purity) ----------------------#
# hin	and hin_jacob are needed for Augmented Lagrangian Minimization 
# Algorithm for optimization
# hin: a vector function specifying inequality constraints such that 
# hin[j] > 0 for all j
# hin.jac: Jacobian of hin

hin_func_noPurity <- function(x, ...){
  return(c(x-0.005, 1-sum(x)-0.005))
}

hin_jacob_noPurity <- function(x, ...){
  return(rbind(diag(1,length(x)), rep(-1,length(x))))
}

logLik_noPurity <- function(x, logY, Z, sigma2C, sigma2A, EM_wgt){
  log_eta_ij = log(Z %*% x)
  mu_ijC = log_eta_ij - sigma2C/2
  mu_ijA = log_eta_ij - sigma2A/2
  
  out = sum(EM_wgt*dnorm(logY, mean = mu_ijC, sd = sqrt(sigma2C), log = TRUE)) + 
    sum((1-EM_wgt)*dnorm(logY, mean = mu_ijA, sd = sqrt(sigma2A), log = TRUE))
  
  return(out)
}

#----------------- Likelihood (given purity) ----------------------#
hin_func_givenPurity <- function(x, rho_i0, ...){
  return(c(x-0.005, 1-rho_i0-sum(x)-0.005))
}

hin_jacob_givenPurity <-function(x, rho_i0, ...){
  return(rbind(diag(1,length(x)), rep(-1,length(x))))
}

logLik_givenPurity <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho  = c(x, 1-rho_i0-sum(x))
  
  eta_ij = Z %*% rho
  mu_ijC = log(eta_ij) - sigma2C/2
  mu_ijA = log(eta_ij) - sigma2A/2
  
  out = sum(EM_wgt*dnorm(logY, mean = mu_ijC, sd = sqrt(sigma2C), log = TRUE)) + 
    sum((1-EM_wgt)*dnorm(logY, mean = mu_ijA, sd = sqrt(sigma2A), log = TRUE))
  
  return(out)
}

#------------------- Actual Update Functions ------------------------#

updatePropn_Single <- function(x, Z, nCT, nG, maxIter_prop, givenPurity, 
                               var_wgt){
  #----------------------------------------#
  # Extract Info                           #
  #----------------------------------------#
  rho_i0 = x[1]
  urho_1 = x[c(2:(nCT+1))]
  logY   = x[c((nCT+2):(nCT+nG+1))]
  EM_wgt = x[c((nCT+nG+2):(nCT+2*nG+1))]
  
  #----------------------------------------#
  # Initialize Values                      #
  #----------------------------------------#
  eta_ij     = drop(Z %*% matrix(urho_1, ncol=1))
  log_eta_ij = log(eta_ij)
  
  sigma2C_1  = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij,
                             EM_wgt = EM_wgt, AB_Up = FALSE, var_wgt=var_wgt)
  sigma2A_1  = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij, 
                            EM_wgt = EM_wgt, AB_Up = TRUE, var_wgt=var_wgt)
  
  if(! is.null(var_wgt)){
    sigma2C_1 = var_wgt*sigma2C_1
    sigma2A_1 = var_wgt*sigma2A_1
  }
  
  mu_ijC = log_eta_ij - sigma2C_1/2
  mu_ijA = log_eta_ij - sigma2A_1/2
  
  logLik_1 = sum(EM_wgt*dnorm(logY,mean=mu_ijC,sd=sqrt(sigma2C_1),log=TRUE))+
    sum((1-EM_wgt)*dnorm(logY,mean=mu_ijA,sd=sqrt(sigma2A_1),log=TRUE))
  
  for(k in 1:maxIter_prop){
    #----------------------------------#
    # Reset Current Values
    #----------------------------------#
    urho_0    = urho_1
    sigma2C_0 = sigma2C_1
    sigma2A_0 = sigma2A_1
    logLik_0  = logLik_1
    
    #----------------------------------#
    # Update Rho values
    #----------------------------------#
    if(givenPurity){
      auglagOut = auglag(par = urho_0[-c(nCT)], fn = logLik_givenPurity, 
                         gr = gradFunc_givenPurity, hin = hin_func_givenPurity, 
                         hin.jac = hin_jacob_givenPurity, logY = logY, 
                         rho_i0 = rho_i0, Z=Z, sigma2C = sigma2C_0, 
                         sigma2A = sigma2A_0, EM_wgt = EM_wgt, 
                         control.optim = list(fnscale=-1), 
                         control.outer = list(trace=FALSE))
      
      urho_1[c(1:(nCT-1))] = auglagOut$par
      urho_1[nCT] = 1-rho_i0-sum(auglagOut$par)
      
      urho_1 = correctRho(est = urho_1, total = 1-rho_i0)
    } else {
      auglagOut = auglag(par = urho_0, fn = logLik_noPurity, 
                         gr = gradFunc_noPurity, hin = hin_func_noPurity, 
                         hin.jac = hin_jacob_noPurity, logY = logY, 
                         Z=Z, sigma2C = sigma2C_0, 
                         sigma2A = sigma2A_0, EM_wgt = EM_wgt, 
                         control.optim=list(fnscale=-1), 
                         control.outer = list(trace=FALSE))
      
      urho_1 = auglagOut$par
      urho_1 = correctRho_v2(urho_1)
    }
    
    #----------------------------------#
    # Update Computational Values
    #----------------------------------#
    eta_ij     = drop(Z %*% matrix(urho_1, ncol=1))
    log_eta_ij = log(eta_ij)
    sigma2C_1s = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij, 
                               EM_wgt = EM_wgt, AB_Up = FALSE, var_wgt=var_wgt)
    sigma2A_1s = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij,
                               EM_wgt = EM_wgt, AB_Up = TRUE, var_wgt=var_wgt)
    
    if(! is.null(var_wgt)){
      sigma2C_1 = var_wgt*sigma2C_1s
      sigma2A_1 = var_wgt*sigma2A_1s
    }else{
      sigma2C_1 = sigma2C_1s
      sigma2A_1 = sigma2A_1s
    }
    
    mu_ijC    = log_eta_ij - sigma2C_1/2
    mu_ijA    = log_eta_ij - sigma2A_1/2
    
    logLik_1  = sum(EM_wgt*dnorm(logY,mean=mu_ijC,sd=sqrt(sigma2C_1),log=TRUE))+
      sum((1-EM_wgt)*dnorm(logY,mean=mu_ijA,sd=sqrt(sigma2A_1),log=TRUE))
    
    #if(logLik_1 < logLik_0 - 0.001){
    #  stop("log likelihood decreased during the EM algorithm.")
    #}
    
    if(max(abs(urho_1-urho_0)) < 1e-3){ break }
  }
  return(list(rho_i=urho_1, iter=k, sigma2C_i=sigma2C_1s, sigma2A_i=sigma2A_1s))
}

#------------------- Update Multiple Subjects ----------------------#
updatePropn_All <- function(logY, rho, tumorPurity, Z, maxIter_prop, 
                            EM_wgt, givenPurity, var_wgt){
  
  PropnInfo = rbind(tumorPurity, rho, logY, EM_wgt)
  
  out = apply(X=PropnInfo, MARGIN = 2, FUN = updatePropn_Single, 
              Z=Z, givenPurity = givenPurity, maxIter_prop = maxIter_prop, 
              nCT=(nrow(rho)), nG = nrow(logY), var_wgt=var_wgt)
  
  #--- Reshape Output ---#
  rho_Mat   = extractFunction(compList = out, element = "rho_i")
  sig2C_Mat = extractFunction(compList = out, element = "sigma2C_i")
  sig2A_Mat = extractFunction(compList = out, element = "sigma2A_i")
  
  return(list(rho_Curr = rho_Mat, sig2C_Curr = sig2C_Mat, 
              sig2A_Curr = sig2A_Mat))
}

#-------------------------------------------------------------------#
#                            EM Weights                             #
#-------------------------------------------------------------------#

updateWgts <- function(logY, rho, sigma2C, sigma2A, Z, propC, var_wgt){
  
  EM_wgt = matrix(NA, nrow=nrow(logY), ncol=ncol(logY))
  
  for(i in 1:ncol(logY)){
    
    logY_i = logY[,i]
    rho_i  = rho[,i]
    
    if(is.null(var_wgt)){
      sigma2C_i = sigma2C[i]
      sigma2A_i = sigma2A[i]
    }else{
      sigma2C_i = sigma2C[i]*var_wgt
      sigma2A_i = sigma2A[i]*var_wgt
    }
    
    # Cmu_ij and Amu_ij are mean values for consistent/aberrant genes
    eta_ij = Z %*% matrix(rho_i, ncol=1)
    Cmu_ij = log(eta_ij) - sigma2C_i/2
    Amu_ij = log(eta_ij) - sigma2A_i/2
    C_lLik = dnorm(logY_i, mean = Cmu_ij, sd = sqrt(sigma2C_i), log = TRUE)
    A_lLik = dnorm(logY_i, mean = Amu_ij, sd = sqrt(sigma2A_i), log = TRUE)
    
    #-----------------------------------#
    # Compiling Weights                 #
    #-----------------------------------#
    EM_wgt[,i] = 1/(1+((1-propC[i])/propC[i])*exp(A_lLik-C_lLik))
  }
  
  return(EM_wgt)
}

PropPlus_Update<- function(Y, rho_0, tumorPurity, givenPurity, 
                           sigma2C_0, sigma2A_0, Z, propC_0,
                           maxIter_PP, maxIter_prop, nG,
                           rhoConverge, var_wgt){
  logY = log(Y)
  
  rho_t1     = rho_0
  propC_t1   = propC_0
  sigma2C_t1 = sigma2C_0
  sigma2A_t1 = sigma2A_0
  
  for(j in 1:maxIter_PP){
    rho_t0     = rho_t1
    propC_t0   = propC_t1
    sigma2C_t0 = sigma2C_t1
    sigma2A_t0 = sigma2A_t1
    
    # Update EM Weights
    EM_wgt = updateWgts(logY = logY, rho = rho_t0, 
                        sigma2C = sigma2C_t0, sigma2A = sigma2A_t0, 
                        Z = Z, propC = propC_t0, var_wgt=var_wgt)
    
    # Update Proportions
    PropP_t1  = updatePropn_All(logY = logY, rho = rho_t0, 
                                tumorPurity = tumorPurity,
                                Z = Z, maxIter_prop = maxIter_prop, 
                                EM_wgt = EM_wgt, givenPurity = givenPurity,
                                var_wgt=var_wgt)
    
    rho_t1     = PropP_t1$rho_Curr
    sigma2C_t1 = PropP_t1$sig2C_Curr
    sigma2A_t1 = PropP_t1$sig2A_Curr
    
    #---- Update proportion of consistent genes
    propC_t1 = colSums(EM_wgt)/nG
    
    #---- Check convergence 
    rho_diff = abs(rho_t0 - rho_t1)
    max_rho_diff = max(c(rho_diff))
    
    if(max_rho_diff < rhoConverge){ break }
    
    if(j %% 10 == 1){
      message("Iter ",j,": max diff of rho: ",max(max_rho_diff))
    }else{
      if(j %% 10 == 0){
        message(".", appendLF = TRUE)
      }else{
        message(".", appendLF = FALSE)
      }
    }
  }
  
  return(list(rho = rho_t1, sigma2C = sigma2C_t1, sigma2A = sigma2A_t1, 
              propC = propC_t1, Iter = j, max_rho_diff = max_rho_diff))
}

#-----------------------------------------------------#
# estimate cell type-specific expression data matrix
#-----------------------------------------------------#

refEstimate <- function(X, borrow4SD=TRUE){

  cellType = colnames(X)
  
  if(any(cellType == "")){
    stop("colnames(X) are cell type labels and cannot be empty.")
  }
  
  if(any(is.na(X))){
    stop("X must not contain any NA entries.")
  }
  
  if(any(X < 0)){
    stop("X should be non-negative.")
  }

  if(any(X < 1e-4)){
    message("Adding 1e-5 to X to ensure valid log transformation.")
    X = X + 1e-5
  }
  
  sortCT  = sort(unique(cellType))
  
  ctTable = table(cellType)
  
  if(min(ctTable) < 2){
    stop("At least two pure samples are needed for each cell type.")
  }
  
  #-----------------------------------------------------#
  # estimate signature matrix                           #
  #-----------------------------------------------------#
  
  logX   = log(X) 
  
  CT_MU  = t(apply(X = logX, MARGIN = 1, FUN = meanFun, 
                   group = cellType))
  
  CT_var = t(apply(X = logX, MARGIN = 1, FUN = varFun,
                   group = cellType))
  
  if(any(colnames(CT_MU) != sortCT) || any(colnames(CT_var) != sortCT)){
    stop("Problems with tapply order!")
  }
  
  refVar = CT_var
  
  # Borrow information across cell types for SD estimation 
  if(borrow4SD){
    X1     = logX
    ntot_p = length(cellType)
    
    for(i in 1:length(sortCT)){
      wwi = which(cellType == sortCT[i])
      X1[,wwi] = logX[,wwi] - CT_MU[,i]
    }
    
    varXPop = apply(X1, 1, var)
    
    for(i in 1:length(sortCT)){
      wi = (sum(cellType == sortCT[i]))/ntot_p
      CT_var[,i] = wi*CT_var[,i] + (1-wi)*varXPop
    }
  }
  
  #-----------------------------------------------------#
  # Generate cell type-specific expressio matrix        #
  #-----------------------------------------------------#
  Z = exp(CT_MU + CT_var/2)
  
  return(list(refMat=Z, refVar=refVar, ctMu=CT_MU, ctVar=CT_var))
}


#-----------------------------------------------------#
# main function                                       #
#-----------------------------------------------------#

ICeDT <- function(Y, Z, tumorPurity=NULL, refVar=NULL, rhoInit=NULL, 
                  maxIter_prop=100, maxIter_PP=100, rhoConverge=1e-3){
  
  #-----------------------------------------------------#
  # Check Input                                         #
  #-----------------------------------------------------#
  
  if(any(duplicated(colnames(Z)))){
    stop("reference matrix Z has duplicated colnames.")
  }
  
  if(any(colnames(Z) == "")){
    stop("colnames(Z) are cell type labels and cannot be empty.")
  }
  
  if(nrow(Y)!=nrow(Z)){
    stop("Y and Z do not have the same number of genes.")
  }

  if(!identical(rownames(Y), rownames(Z))){
    stop("rownames of Y and Z are different.")
  }
  
  if(!is.null(refVar)){
    if(!identical(dim(Z), dim(refVar))){
      stop("the dimensions of Z and refVar do not match.")
    }
    
    if(!identical(rownames(Z), rownames(refVar))){
      stop("row names of Z and refVar do not match.")
    }
    
    if(!identical(colnames(Z), colnames(refVar))){
      stop("column names of Z and refVar do not match.")
    }
  }
  
  if(any(is.na(Y)) | any(is.na(Z))){
     stop("Y and Z must not contain any NA entries.")
  }
  
  if(any(Y < 0) | any(Z < 0)){
    stop("Y and Z should be non-negative.")
  }
  
  if(any(Y < 1e-4)){
    message("Adding 1e-5 to Y to ensure valid log transformation.")
    Y = Y + 1e-5
  }
  
  if(any(Z < 1e-4)){
    message("Adding 1e-5 to Z to ensure valid log transformation.")
    Z = Z + 1e-5
  }
  
  nG  = nrow(Y) # number of genes
  nS  = ncol(Y) # number of bulk tumor (mixed cell type) samples
  nCT = ncol(Z) # number of cell types
  
  if(! is.null(rhoInit)){
    if(nrow(rhoInit) != nCT){
      stop("nrow(rhoInit) should equal to the number of cell types.")
    }
    
    if(ncol(rhoInit) != nS){
      stop("ncol(rhoInit) should equal to the number bulk tumor samples.")
    }
    
    if(! setequal(rownames(rhoInit), colnames(Z))){
      stop("rownames of rhoInit (cell types) do not match with Z.")
    }else{
      rhoInit = rhoInit[match(colnames(Z),rownames(rhoInit)),]
    }
    
    if(any(colnames(rhoInit) != colnames(Y))){
      stop("colnames of rhoInit (samples) do not match with Y.")
    }
  }
  
  givenPurity = !(is.null(tumorPurity))

  if(givenPurity){
    if(ncol(Y) != length(tumorPurity)){
      stop("ncol of Y does not equal to the length of tumorPurity.")
    } else if(!all(colnames(Y) == names(tumorPurity))){
      stop("colnames of Y do not match with the names of tumorPurity.")
    }
  }else{
    tumorPurity = rep(0, ncol(Y))
    names(tumorPurity) = colnames(Y)
  }
  
  #-----------------------------------------------------#
  # initialization                                      #
  #-----------------------------------------------------#
  
  if(! is.null(refVar)){
    maxVar  = apply(X = refVar, MARGIN = 1, FUN = max)
    var_wgt = maxVar/median(maxVar)
    
    quants = quantile(var_wgt,prob=c(0.15, 0.85),na.rm=TRUE)
    wlo = which(var_wgt < quants[1])
    wup = which(var_wgt > quants[2])
    var_wgt[wlo] = quants[1]
    var_wgt[wup] = quants[2]
  }else{
    var_wgt = NULL
  }
  
  # percent of consistent genes per sample 
  propC_0 = rep(0.5, nS)
  
  # lmInit function estimate cell type compositon by linear regression
  if(is.null(rhoInit)){
    rho_0 = apply(X = rbind(tumorPurity,Y), MARGIN=2, FUN=lmInit, Z=Z)
  }else{
    if(givenPurity){
      rho_0 = t(t(rhoInit)/colSums(rhoInit) * (1 - tumorPurity))
    }else{
      rho_0 = rhoInit
    }
  }

  # estimate residual variance for consistent/aberrant genes
  sigma_0 = apply(X = rbind(Y, rho_0), MARGIN=2, FUN=sigmaInit, Z=Z, 
                  nG=nG, var_wgt=var_wgt)
  
  sigma2C_0 = sigma_0[1,]
  sigma2A_0 = sigma_0[2,]
  
  #-----------------------------------------------------#
  # EM algorithm                                        #
  #-----------------------------------------------------#
  PropPlus_Out = suppressWarnings(
    PropPlus_Update(Y = Y, rho_0 = rho_0, tumorPurity = tumorPurity, 
                    givenPurity = givenPurity, sigma2C_0 = sigma2C_0, 
                    sigma2A_0 = sigma2A_0, Z = Z, propC_0 = propC_0, 
                    maxIter_PP = maxIter_PP, maxIter_prop = maxIter_prop, 
                    nG = nG, rhoConverge = rhoConverge, var_wgt = var_wgt))
  
  rho_1     = PropPlus_Out$rho
  sigma2C_1 = PropPlus_Out$sigma2C
  sigma2A_1 = PropPlus_Out$sigma2A
  propC_1   = PropPlus_Out$propC
  
  EM_wgt =  updateWgts(logY = log(Y), rho = rho_1, 
                       sigma2C = sigma2C_1, sigma2A = sigma2A_1, 
                       Z = Z, propC = propC_1, var_wgt=var_wgt)
  
  #-----------------------------------------------------#
  # OUTPUT                                              #
  #-----------------------------------------------------#
  
  if(givenPurity){
    rho_final  = rbind(tumorPurity, rho_1)
  } else {
    tumorPurity_est = 1 - colSums(rho_1)
    rho_final  = rbind(tumorPurity_est, rho_1)
  }
  
  rownames(rho_final) = c("tumor", colnames(Z))
  colnames(rho_final) = colnames(Y)
  
  w2switch = which(sigma2C_1 > sigma2A_1)
  if(length(w2switch) > 0){
    sigma_tmp           = sigma2A_1[w2switch]
    sigma2A_1[w2switch] = sigma2C_1[w2switch]
    sigma2C_1[w2switch] = sigma_tmp
    propC_1[w2switch]   = 1 - propC_1[w2switch]
    EM_wgt[,w2switch]   = 1 - EM_wgt[,w2switch]
  }
  
  outList = list(rho     = rho_final,
                 sigma2C = sigma2C_1,
                 sigma2A = sigma2A_1,
                 cProp   = propC_1,
                 cProb   = EM_wgt)
  
  return(outList)
}

