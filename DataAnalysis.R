### Necessary packages
library(survival)
library(glmnet)
library(grpreg)
library(gglasso)
library(haven)

### Loading Data
GBSG2 <- read_dta("rott2.dta")
View(GBSG2)
str(GBSG2)

### Different Analyses (without interactions)
# Cox model without confounders (unadjusted analysis)
coxph(Surv(rf, rfi) ~ chemo, 
      GBSG2, robust=TRUE)
# 0.17347 (0.06092)

# Cox model with all main effects
coxph(Surv(rf, rfi) ~ chemo+year+age+meno+as.factor(size)+grade+nodes+pr+er+hormon, 
      GBSG2, robust=TRUE)
# -0.1239273  (0.0749121)

# Lasso
modelMatrix=model.matrix(rf~ -1 + chemo+year+age+meno+as.factor(size)+grade+nodes+pr+er, data=GBSG2)
modelMatrix=modelMatrix[,-5]

set.seed(123)
group=factor(c("chemo", "year", "age", "meno", "size", "size", "grade", "nodes", "pr", "er"))
cv.lasso=cv.grpsurv(modelMatrix, Surv(GBSG2$rf, GBSG2$rfi), group, nfolds=20, alpha = 1)
ind=which(cv.lasso$lambda==cv.lasso$lambda.min)
ind2=min(which(cv.lasso$cve<=(cv.lasso$cve[ind]+cv.lasso$cvse[ind])))

model_cox = grpsurv(modelMatrix, Surv(GBSG2$rf, GBSG2$rfi), group, lambda=c(0.05, cv.lasso$lambda[ind2]),
                    alpha = 1)

selection1 <- which(c(1,model_cox$beta[2:10, 2] != 0) != 0) ## always includes the treatment


model1.cox = coxph(Surv(GBSG2$rf, GBSG2$rfi) ~ modelMatrix[,selection1], robust=TRUE)
summary(model1.cox)
#-0.020897 (0.066240)
# chemo, size (factor), grade and nodes
alpha = coef(model1.cox)[1]
betas = coef(model1.cox)[-1]

# Poor man's approach
cv.lasso=cv.grpsurv(modelMatrix, Surv(GBSG2$rf, 1-GBSG2$rfi), group, nfolds=20, alpha = 1)
ind=which(cv.lasso$lambda==cv.lasso$lambda.min)
ind2=min(which(cv.lasso$cve<=(cv.lasso$cve[ind]+cv.lasso$cvse[ind])))

model = grpsurv(modelMatrix, Surv(GBSG2$rf, 1-GBSG2$rfi), group, lambda=c(0.05, cv.lasso$lambda[ind2]),
                alpha = 1)
selection2 = which( model$beta[1:10,2] != 0)
# chemo, year and age


trt = cv.grpreg(modelMatrix[,2:10], modelMatrix[,1], group[-1],
                nfolds=20, family="binomial", intercept=TRUE, standardize=FALSE, alpha = 1)


ind=which(trt$lambda==trt$lambda.min)
ind2=min(which(trt$cve<=(trt$cve[ind]+trt$cvse[ind])))

model = grpreg(modelMatrix[,2:10], modelMatrix[,1], group[-1], lambda=c(0.05, trt$lambda[ind2]),
               alpha = 1,intercept=TRUE, family="binomial")

selection_trt = which(model$beta[-1,2] != 0)+1 
# age, meno, nodes

selection_PMA = union(selection1, union(selection2, selection_trt))
model_PMA = coxph(Surv(GBSG2$rf, GBSG2$rfi) ~ modelMatrix[,selection_PMA], robust=TRUE)
summary(model_PMA)
#-0.127579 (0.074421)


# Triple Selection Approach
la  = rep(0,2982)                           # Gradient w.r.t parameter of interest
lb  = matrix(rep(0,9*2982),nrow = 2982) 

stime = sort(GBSG2$rf)          # Sorted survival/censored times
otime = order(GBSG2$rf)         # Order of time

status=GBSG2$rfi
bhaz = c()
sS0 = c()
sS1 = matrix(rep(0,(9+1)*2982),nrow = 2982)
k   = 1
while( k<=2982)
{
  if (status[otime[k]]==1)
  {
    
    ind = which(GBSG2$rf >= stime[k]) 
    S0  = 0
    S1  = rep(0,9+1)
    
    if (length(ind)>0)
    {
      for (j in 1:length(ind))
      {
        tmp = exp(modelMatrix[,1][ind[j]]*alpha+
                    modelMatrix[,2:10][ind[j],selection1[-1]-1]%*%betas)
        S0  = S0 + tmp
        
        S1 = S1 + tmp %*%t(modelMatrix[ind[j],])
        
      }
    }
    S0 = apply(S0,1,as.numeric)
    
    bhaz[k] = sum(as.numeric(GBSG2$rf==stime[k]&status==1))/S0
    sS0[k] = S0
    sS1[k,] = S1
    la[k]  = -(modelMatrix[,1][otime[k]] - S1[1]/S0)
    lb[k,] = -(modelMatrix[,2:10][otime[k],c(1:9)] - S1[c(2:(9+1))]/S0)
  }
  k = k + 1
  
}

group_num = c(1, 2, 3, 4, 5, 5, 6, 7, 8, 9)

fit.cv <- cv.gglasso(lb, la, group_num[-1]-1, intercept = FALSE, nfolds=20)
best_lambda_fit.cv <- fit.cv$lambda.1se
fit <- gglasso(lb, la, group_num[-1]-1, loss = "ls", intercept = FALSE, lambda = best_lambda_fit.cv)
selection3 = which(fit$beta != 0)+1 
# age, nodes, pr and er 

selection_prop = union(selection1, union(selection2, selection3))
model_prop = coxph(Surv(GBSG2$rf, GBSG2$rfi) ~ modelMatrix[,selection_prop], robust=TRUE)
summary(model_prop)
#-1.441e-01 (7.398e-02)
