### Necessary Packages

library("glmnet")
library("survival")
library("mvtnorm")
library("glasso")
library("lpSolve")
library("MASS")

### Code from the Fang et al (2017) paper
coxhdi2 <- function(x,time,status,method1="S",cptw="L",indx = c(1))
{
  ############################################################################################################# 
  ##Input values:                                                                                            ##
  ##x:       Covariates                                                                                       ##
  ##time:    Survival times or censoring times                                                             ##
  ##status:  Censoring status
  ##method1: "S","W" or "L", represents Score, Wald or Likelihood Ratio Test                                 ##
  ##cptw:    "L" or "D", represents using either Lasso or Dantzig method to compute w (decorrelation vector),##
  ##         where Lasso gives faster computation                                                            ##
  ##indx:    Set of indices to conduct testing
  #############################################################################################################
  
  Pval = c()
  d    = dim(x)[2]
  n    = dim(x)[1]
  for (coi in indx)
  {
    lambdas = seq(from = 0.1,to=1.5,by = 0.03) 
    lambdas = lambdas*log(d)/n                 ##Set of tuning parameters
    cv.fit <- cv.glmnet(x, Surv(time, status),family = "cox", lambda = lambdas, maxit = 10000,nfolds = 5,alpha=1)
    fit     = cv.fit$glmnet.fit
    tmp     = which(fit$lambda == cv.fit$lambda.min)
    
    if (sum(fit$beta[,tmp]!=0)==0)
    {
      beta = rep(0,d)
    } else {
      beta    = fit$beta[,tmp]   # Initial Estimator
    }
    betas   = beta
    #betas[coi] = 0              # H0: beta_i = 0
    
    stime = sort(time)          # Sorted survival/censored times
    otime = order(time)         # Order of time
    
    Vs  = matrix(rep(0,d*d),nrow = d)
    Hs  = Vs                                 # Hessian
    ind = 0
    
    la  = rep(0,n)                           # Gradient w.r.t parameter of interest
    lb  = matrix(rep(0,(d-1)*n),nrow = n)    # Gradient w.r.t nuisance parameter (theta)
    i   = 1
    while( i<=n)
    {
      if (status[otime[i]]==1)
      {
        ind = which(time >= stime[i])
        S0  = 0
        S1  = rep(0,d)
        S2  = matrix(rep(0,d*d),nrow = d)
        
        if (length(ind)>0)
        {
          for (j in 1:length(ind))
          {
            tmp = exp(x[ind[j],]%*%betas)
            S0  = S0 + tmp
            
            S1  = S1 + tmp %*%t(x[ind[j],])
            
            tmp = apply(tmp,1,as.numeric)   	
            S2  = S2 + tmp*x[ind[j],]%*%t(x[ind[j],])
          }
        }
        S0 = apply(S0,1,as.numeric)
        
        la[i]  = -(x[otime[i],coi] - S1[coi]/S0)
        if (coi == 1)
        {
          lb[i,] = -(x[otime[i],c((coi+1):d)] - S1[c((coi+1):d)]/S0)
        } else if (coi == d){
          lb[i,] = -(x[otime[i],c(1:(coi-1))] - S1[c(1:(coi-1))]/S0)
        } else {
          lb[i,] = -(x[otime[i],c(1:(coi-1), (coi+1):d)] - S1[c(1:(coi-1), (coi+1):d)]/S0)
        }
        V   = S0*S2 - t(S1)%*%(S1)
        Hs  = Hs + V/(n*S0^2)          
      }
      i = i + 1
    }
    
    #Hs = Hs/n
    
    if (cptw == "L")
    {
      if (method1=="L")
      {
        cv.fit = glmnet(lb,la,standardize = FALSE,intercept = FALSE)
        fit     = cv.fit$glmnet.fit
        tmp     = which(fit$lambda == cv.fit$lambda.min)
        if (sum(fit$beta[,tmp]!=0)==0)
        {
          what = rep(0,d-1)
        } else {
          what    = fit$beta[,tmp]   
        }
      } else {
        fit    = glmnet(lb,la,standardize = FALSE,intercept = FALSE)
        lgth = length(fit$lambda)
        HBIC = rep(10000,lgth)
        for (tmp in 1:lgth)
        {   
          if (sum(fit$beta[,tmp])!=0)
          {
            loss = sum((la-lb%*%fit$beta[,tmp])^2)/n
            HBIC[tmp] = log(loss) + fit$df[tmp]*log(d)/n*eta*(log(log(n)))
          }
        }
        tmp = which(HBIC==min(HBIC))
        what = fit$beta[,tmp]
      }
      
      
    } else if (cptw == "D") {
      if (coi == 1)
      {
        Hab = Hs[(coi+1):d,coi]
        Hbb = Hs[(coi+1):d,(coi+1):d]
      } else if (coi == d){
        Hab = Hs[1:(coi-1),coi]
        Hbb = Hs[1:(coi-1),1:(coi-1)]
      } else{
        Hab = Hs[c(1:(coi-1),(coi+1):d),coi]
        Hbb = Hs[c(1:(coi-1),(coi+1):d),c(1:(coi-1),(coi+1):d)]
      }
      
      tmp = matrix(rep(0,(d-1)*(d-1)),nrow = d-1)
      A1  = cbind(Hbb,tmp)
      A2  = cbind(-Hbb,tmp)
      A3  = cbind(diag(d-1),-diag(d-1))
      A4  = cbind(-diag(d-1),-diag(d-1))
      
      A   = rbind(A1,A2,A3,A4)
      obj = c(rep(0,d-1),rep(1,d-1))
      dir = rep("<=",(d-1)*2 + (d-1)*2)
      lambda = 1*sqrt(log(d)/n)
      rhs = c(Hab + lambda, -Hab + lambda, rep(0,d-1),rep(0,d-1))
      tmp = lp(direction = "min",objective.in = obj,const.mat = A,const.dir = dir, const.rhs = rhs)
      what = tmp$solution[1:(d-1)]
    }
    
    what = matrix(what,nrow=d-1)
    if (method1 == "S"){
      # Decorrelated Score
      if (coi == 1)
      {
        var   = max(Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi],0.1)
      } else if (coi == d){
        var   = max(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi],0.1)
      } else {
        var   = max(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi],0.1)
      }
      lbt = colSums(lb)
      lbt = matrix(lbt,nrow = d-1)
      S = sum(la) - t(what)%*%(lbt)
      S = S/n
      
      tmp   = sqrt(n)*S/sqrt(var)
      pval  = 2*pnorm(-abs(tmp))
    } else if (method1 == "W") {
      # Decorrelated Wald
      if (coi == 1)
      {
        S = beta[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi])
        var   = Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi]
      } else if (coi == d){
        S = beta[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi])
        var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi]
      } else {
        S = beta[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi])
        var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi]
      }
      tmp   = (n)*S^2*(max(var,1e-8))
      pval  = 1-pchisq(tmp,1)
    } else if (method1 == "L") {
      # Deccorelated LHR
      tmp = which(la!=0)
      
      
      if (coi == 1)
      {
        atilde = beta[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi])
        beta1 = c(atilde,beta[(coi+1):d]-atilde*what[(coi):(d-1)])
      } else if (coi == d) {
        atilde = beta[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi])
        beta1 = c(beta[1:(coi-1)]-atilde*what[1:(coi-1)],atilde)
      } else {
        atilde = beta[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi])
        beta1 = c(beta[1:(coi-1)]-atilde*what[1:(coi-1)],atilde,beta[(coi+1):d]-atilde*what[(coi):(d-1)])
      }
      
      
      lossa  = 0    #Log-Likelihood under alternative
      i = 1
      while( i<=n)
      {
        ind = which(time >= stime[i])
        S0  = 0
        # S1  = matrix(rep(0,d),nrow = d)
        if (length(ind)>0)
        {
          for (j in 1:length(ind))
          {
            tmp = 1/n * exp(x[ind[j],]%*%beta1)
            tmp = apply(tmp,1,as.numeric)
            S0  = S0 + tmp
          }
        }
        if (status[otime[i]]==1)
        {
          lossa = lossa + 1/n*(x[otime[i],]%*%beta1 - log(S0))
        }
        i = i + 1
      }
      lossn  = 0  #Log-Likelihood under null
      i = 1
      if (coi == 1)
      {
        beta2 = c(0,beta[(coi+1):d])
      } else if (coi == d) {
        beta2 = c(beta[1:(coi-1)],0)
      } else {
        beta2 = c(beta[1:(coi-1)],0,beta[(coi+1):d])
      }
      while( i<=n)
      {
        ind = which(time >= stime[i])
        S0  = 0
        if (length(ind)>0)
        {
          for (j in 1:length(ind))
          {
            tmp = 1/n * exp(x[ind[j],]%*%beta2)
            tmp = apply(tmp,1,as.numeric) 
            S0  = S0 + tmp
          }
        }
        
        if (status[otime[i]]==1)
        {
          lossn = lossn + 1/n*(x[otime[i],]%*%beta2 - log(S0))
        }
        i = i + 1
      }
      tmp   = n*lossa - n*lossn
      tmp   = 2*tmp
      pval  = 1-pchisq(tmp,1)
    }
    Pval = c(Pval,pval)
  }
  return(Pval)
}

coxhdiHat <- function(x,time,status,fit_cox)
{
  ############################################################################################################# 
  ##Input values:                                                                                            ##
  ##x:       Covariates                                                                                       ##
  ##time:    Survival times or censoring times                                                             ##
  ##status:  Censoring status
  ##method1: "S","W" or "L", represents Score, Wald or Likelihood Ratio Test                                 ##
  ##cptw:    "L" or "D", represents using either Lasso or Dantzig method to compute w (decorrelation vector),##
  ##         where Lasso gives faster computation                                                            ##
  ##indx:    Set of indices to conduct testing
  #############################################################################################################
  
  d    = dim(x)[2]
  n    = dim(x)[1]
  
  
  #cv.fit = cv.glmnet(x, Surv(time, status),family = "cox", maxit = 10000, nfolds = 20,alpha=1)
  #fit = glmnet(x, Surv(time, status),family = "cox",alpha=1, lambda = cv.fit$lambda.1se)
  
  fit=fit_cox
  
  if (sum(fit$beta!=0)==0)
  {
    beta = rep(0,d)
  } else {
    beta    = fit$beta  # Initial Estimator
  }
  betas   = beta
  #betas[1] = 0              # H0: beta_i = 0
  
  stime = sort(time)          # Sorted survival/censored times
  otime = order(time)         # Order of time
  
  Vs  = matrix(rep(0,d*d),nrow = d)
  Hs  = Vs                                 # Hessian
  ind = 0
  
  la  = rep(0,n)                           # Gradient w.r.t parameter of interest
  lb  = matrix(rep(0,(d-1)*n),nrow = n)    # Gradient w.r.t nuisance parameter (theta)
  
  select=which(betas!=0)
  betas_select=betas[select]
  
  
  i   = 1
  while( i<=n)
  {
    if (status[otime[i]]==1)
    {
      ind = which(time >= stime[i])
      S0  = 0
      S1  = rep(0,d)
      S2  = matrix(rep(0,d*d),nrow = d)
      
      if (length(ind)>0)
      {
        for (j in 1:length(ind))
        {
          tmp = exp(x[ind[j],select]%*%betas_select)
          S0  = S0 + tmp
          
          S1  = S1 + tmp %*%t(x[ind[j],])
          
          tmp = apply(tmp,1,as.numeric)   	
          S2  = S2 + tmp*x[ind[j],]%*%t(x[ind[j],])
        }
      }
      S0 = apply(S0,1,as.numeric)
      
      la[i]  = -(x[otime[i],1] - S1[1]/S0)
      
      lb[i,] = -(x[otime[i],c((1+1):d)] - S1[c((1+1):d)]/S0)
      
      V   = S0*S2 - t(S1)%*%(S1)
      Hs  = Hs + V/(n*S0^2)          
    }
    i = i + 1
  }
  
  #Hs = Hs/n
  
  
  cv.lasso = cv.glmnet(lb, la, alpha = 1, nfolds=20, intercept=FALSE, standardize=FALSE)
  fit  = glmnet(lb, la, alpha = 1, lambda = cv.lasso$lambda.1se, intercept=FALSE, standardize=FALSE)
  
  what = fit$beta
  what = matrix(what,nrow=d-1)
  
  # Decorrelated Wald
  
  S = beta[1] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[1,1] - t(what)%*%Hs[c((1+1):d),1])
  var   = Hs[1,1] - t(what)%*%Hs[c((1+1):d),1]
  
  tmp   = (n)*S^2*(max(var,1e-8))
  pval  = 1-pchisq(tmp,1)
  
  return(pval)
  
}

### Function that simulates data and evaluate the data using different
# estimation and variable selection approaches
simFunction=function(c_beta, gamma1, c_gamma2){
  
  n=200 # Change sample size for different settings
  p=250 # Change number of variables for different settings
  
  pVal_prop=c(); pVal_PMA=c(); pVal_fangDS=c();
  pVal_DS=c(); pVal_DS2=c(); pVal_LRT=c(); pVal_lasso=c();
  pVal_Fang=c(); pVal_FangHat=c()
  
  
  NumberSelected=c(); TotalSelection=c()
  
  sigma=toeplitz(0.5^(0:p))
  mu = rep(0, p+1)
  
  
  for(i in 1:1000){
    
    set.seed(i)
    
    # Code below simulates the data
    # User can use this code to change data-generating mechanism
    X=mvrnorm(n, mu, sigma)
    D = X[,1]
    X=X[,-1]
    
    p2=10
    beta=c_beta*c( (1:(p2))^(-1), rep(0,p-p2) ) 
    Time = rexp(n, exp(0*D+X%*%beta))
    
    p2=5
    gamma2=c_gamma2*c(  (1:p2)^(-1), (1:p2)^(-1), rep(0,p-2*p2))
    TimeC = rexp(n, exp(gamma1*D+X%*%gamma2))
    
    
    C = TimeC <= Time
    Time[which(C)] = TimeC[which(C)]
    Y = Surv(Time, 1-C)
    
    ####################################################
    
    ### Ordinary method, witch penelization
    cv.lasso = cv.glmnet(cbind(D,X), Y, alpha = 1, family = "cox", nfolds=20)
    model_cox = glmnet(cbind(D,X), Y, alpha = 1, family = "cox", lambda = cv.lasso$lambda.1se)
    selection1 <- which(c(1,model_cox$beta[2:(p+1)] != 0) != 0) ## always includes the treatment
    
    model1.cox = coxph(Y ~ cbind(D,X)[,selection1], robust=TRUE)
    alpha = coef(model1.cox)[1]
    betas = coef(model1.cox)[-1]
    
    
    ### Cox regression on censored events
    Y2 = Surv(Time, C)
    cv.lasso = cv.glmnet(cbind(D,X), Y2, alpha = 1, family = "cox", nfolds=20)
    model = glmnet(cbind(D,X), Y2, alpha = 1, family = "cox", lambda = cv.lasso$lambda.1se)
    selection2 = which( model$beta[1:(p+1)] != 0)
    
    ### Exposure model (biased-reduced method)
    la  = rep(0,n)                           # Gradient w.r.t parameter of interest
    lb  = matrix(rep(0,(p)*n),nrow = n) 
    
    stime = sort(Time)          # Sorted survival/censored times
    otime = order(Time)         # Order of time
    
    status=1-C
    bhaz = c()
    sS0 = c()
    sS1 = matrix(rep(0,(p+1)*n),nrow = n)
    k   = 1
    while( k<=n)
    {
      if (status[otime[k]]==1)
      {
        
        ind = which(Time >= stime[k]) 
        S0  = 0
        S1  = rep(0,p+1)
        
        if (length(ind)>0)
        {
          for (j in 1:length(ind))
          {
            tmp = exp(D[ind[j]]*alpha+X[ind[j],selection1[-1]-1]%*%betas)
            S0  = S0 + tmp
            
            S1 = S1 + tmp %*%t(cbind(D,X)[ind[j],])
            
          }
        }
        S0 = apply(S0,1,as.numeric)
        
        bhaz[k] = sum(as.numeric(Time==stime[k]&status==1))/S0
        sS0[k] = S0
        sS1[k,] = S1
        la[k]  = -(D[otime[k]] - S1[1]/S0)
        lb[k,] = -(X[otime[k],c(1:p)] - S1[c(2:(p+1))]/S0)
      }
      k = k + 1
      
    }
    
    
    cv.lasso = cv.glmnet(lb, la, alpha = 1, nfolds=20, intercept=FALSE, standardize=FALSE)
    model = glmnet(lb, la, alpha = 1, lambda = cv.lasso$lambda.1se, intercept=FALSE, standardize=FALSE)
    selection3 = which(model$beta != 0)+1  
    
    
    ### Exposure Model (poor man's approach)
    trt = cv.glmnet(X, D, alpha = 1, nfolds=20, intercept=TRUE, standardize=TRUE)
    model = glmnet(X, D, alpha = 1, lambda = trt$lambda.1se, intercept=TRUE, standardize=TRUE)
    selection_trt = which(model$beta != 0)
    
    
    ### Final Models
    #Proposal
    selection_prop = union(selection1, union(selection2, selection3))
    model_prop = coxph(Y ~ cbind(D,X)[,selection_prop], robust=TRUE)
    
    #Poor man's approach
    selection_PMA = union(selection1, union(selection2, selection_trt))
    model_PMA = coxph(Y ~ cbind(D,X)[,selection_PMA], robust=TRUE)
    
    #Double Selection based on Fang's method
    selection_fangDS = union(selection1, selection3)
    model_fangDS = coxph(Y ~ cbind(D,X)[,selection_fangDS], robust=TRUE)
    
    #Double Selection (without exposure model)
    selection_DS = union(selection1, selection2)
    model_DS = coxph(Y ~ cbind(D,X)[,selection_DS], robust=TRUE) 
    
    #Double Selection (without censoring model)
    selection_DS2 = union(selection1, selection_trt)
    model_DS2 = coxph(Y ~ cbind(D,X)[,selection_DS2], robust=TRUE)
    
    
    #Logrank test
    model_LRT = coxph(Y ~ D, robust=TRUE)
    
    
    ### Results
    pVal_prop[i]=as.double(summary(model_prop)$coefficients[1,6])
    pVal_PMA[i]=as.double(summary(model_PMA)$coefficients[1,6])
    pVal_fangDS[i]=as.double(summary(model_fangDS)$coefficients[1,6])
    pVal_DS[i]=as.double(summary(model_DS)$coefficients[1,6])
    pVal_DS2[i]=as.double(summary(model_DS2)$coefficients[1,6])
    pVal_LRT[i]=as.double(summary(model_LRT)$coefficients[1,6])
    
    pVal_lasso[i]=as.double(summary(model1.cox)$coefficients[1,6])
    
    pVal_Fang[i] = coxhdi2(cbind(D,X),Time,1-C,method1="W",cptw="D",indx = c(1))
    pVal_FangHat[i] = coxhdiHat(cbind(D,X),Time,1-C, model_cox)
    
  }
  
  list(cbind(pVal_prop, pVal_PMA, pVal_fangDS,
             pVal_DS, pVal_DS2, pVal_LRT, pVal_lasso,
             pVal_Fang, pVal_FangHat))
  
}



### Example
simFunction(c_beta = 0, gamma1 = 1, c_gamma2=0)
