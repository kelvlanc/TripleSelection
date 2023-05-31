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
# Post-Lasso with interactions and higher-order terms
modelMatrix=model.matrix(rf~ -1 + chemo+(as.factor(size)+grade+nodes)^2+I(grade^2)+I(nodes^2), data=GBSG2)
modelMatrix=modelMatrix[,-2]

set.seed(123)
group=factor(c("chemo", "0", "0", "0", "0", "grade2", "nodes2", "sizeGrade", "sizeGrade", "sizeNodes", "sizeNodes", "gradeNodes"))
cv.lasso=cv.grpsurv(modelMatrix, Surv(GBSG2$rf, GBSG2$rfi), group, nfolds=20, alpha = 1)
ind=which(cv.lasso$lambda==cv.lasso$lambda.min)
ind2=min(which(cv.lasso$cve<=(cv.lasso$cve[ind]+cv.lasso$cvse[ind])))

model_cox = grpsurv(modelMatrix, Surv(GBSG2$rf, GBSG2$rfi), group, lambda=c(0.05, cv.lasso$lambda[ind2]),
                    alpha = 1)

selection1 <- which(c(1,model_cox$beta[2:12, 2] != 0) != 0) ## always includes the treatment
modelMatrix_outcome = modelMatrix[,selection1]

model1.cox = coxph(Surv(GBSG2$rf, GBSG2$rfi) ~ modelMatrix[,selection1], robust=TRUE)
summary(model1.cox)
# chemo, size (factor), grade, nodes and nodes^2
# -0.1119695 (0.0646724)
alpha = coef(model1.cox)[1]
betas = coef(model1.cox)[-1]

# Poor man's approach
modelMatrix=model.matrix(rf~ -1 + chemo+(year+age)^2+I(year^2)+I(age^2), data=GBSG2)
set.seed(123)
cv.lasso = cv.glmnet(modelMatrix, Surv(GBSG2$rf, 1-GBSG2$rfi), alpha = 1, family = "cox", nfolds=20, penalty.factor=c(1, rep(0,2), rep(1,3)))
model_cox = glmnet(modelMatrix, Surv(GBSG2$rf, 1-GBSG2$rfi), alpha = 1, family = "cox", lambda = cv.lasso$lambda.1se, penalty.factor=c(1, rep(0,2), rep(1,3)))
selection2 <- which(c(1,model_cox$beta[2:6] != 0) != 0) 
# year and age

modelMatrix=model.matrix(rf~ -1 + chemo+(age+meno+nodes)^2+I(age^2)+I(nodes^2), data=GBSG2)
cv.lasso = cv.glmnet(modelMatrix[,2:9], modelMatrix[,1], alpha = 1, nfolds=20, intercept=TRUE, standardize=FALSE, penalty.factor=c(rep(0,3), rep(1,5)))
fit  = glmnet(modelMatrix[,2:9], modelMatrix[,1], alpha = 1, lambda = cv.lasso$lambda.1se, intercept=TRUE, standardize=FALSE, penalty.factor=c(rep(0,3), rep(1,5)))
fit$beta
# age, meno, nodes, nodes2, age*nodes

modelMatrix_PMA=model.matrix(rf~ -1 + chemo+as.factor(size)+grade+nodes+I(nodes^2)+year+age+meno+age*nodes, data=GBSG2)
modelMatrix_PMA=modelMatrix_PMA[,-2]

model_PMA = coxph(Surv(GBSG2$rf, GBSG2$rfi) ~ modelMatrix_PMA, robust=TRUE)
summary(model_PMA)
#-2.908e-01 (7.717e-02)

# Triple Selection Approach
set.seed(123)
modelMatrix=model.matrix(rf~ -1 +chemo + (age+nodes+ pr+ er)^2 + I(age^2) + I(nodes^2) + I(pr^2) + I(er^2), data=GBSG2)

la  = rep(0,2982)                           # Gradient w.r.t parameter of interest
lb  = matrix(rep(0,14*2982),nrow = 2982) 

stime = sort(GBSG2$rf)          # Sorted survival/censored times
otime = order(GBSG2$rf)         # Order of time

status=GBSG2$rfi
bhaz = c()
sS0 = c()
sS1 = matrix(rep(0,(14+1)*2982),nrow = 2982)
k   = 1
while( k<=2982)
{
  if (status[otime[k]]==1)
  {
    
    ind = which(GBSG2$rf >= stime[k]) 
    S0  = 0
    S1  = rep(0,14+1)
    
    if (length(ind)>0)
    {
      for (j in 1:length(ind))
      {
        tmp = exp(modelMatrix_outcome[,1][ind[j]]*alpha+
                    modelMatrix_outcome[,2:6][ind[j],]%*%betas)
        S0  = S0 + tmp
        
        S1 = S1 + tmp %*%t(modelMatrix[ind[j],])
        
      }
    }
    S0 = apply(S0,1,as.numeric)
    
    bhaz[k] = sum(as.numeric(GBSG2$rf==stime[k]&status==1))/S0
    sS0[k] = S0
    sS1[k,] = S1
    la[k]  = -(modelMatrix[,1][otime[k]] - S1[1]/S0)
    lb[k,] = -(modelMatrix[,2:15][otime[k],c(1:14)] - S1[c(2:(14+1))]/S0)
  }
  k = k + 1
  
}

head(modelMatrix)
cv.lasso = cv.glmnet(lb, la, alpha = 1, nfolds=20, intercept=FALSE, standardize=FALSE, penalty.factor=c(rep(0,4), rep(1,10)))
fit  = glmnet(lb, la, alpha = 1, lambda = cv.lasso$lambda.1se, intercept=FALSE, standardize=FALSE, penalty.factor=c(rep(0,4), rep(1,10)))

fit$beta
# age, nodes, pr, er and pr2

modelMatrix_TS=model.matrix(rf~ -1 + chemo+as.factor(size)+grade+nodes+I(nodes^2)+year+age + pr+er+I(pr^2), data=GBSG2)
modelMatrix_TS=modelMatrix_TS[,-2]

model_TS = coxph(Surv(GBSG2$rf, GBSG2$rfi) ~ modelMatrix_TS, robust=TRUE)
summary(model_TS)
# -3.038e-01 (7.430e-02)