########################################################################################################
## Script that generates synthetic data and performs the best model selection procedure 
## Author: Pilar González-Barquero
#######################################################################################################

####################################
# 1 - Load all needed libraries
####################################
library(glmnet)
library(survival)
library(randomForestSRC)
library(stats)
library(mvtnorm)
library(expint)
library(survAUC)

####################################
# 2 - Load all needed functions
####################################

source("functions/generate_uncensored_data.R")
source("functions/simulate_scale_censoring_proportion.R")
source("functions/generate_censored_data.R")
source("functions/adaptive_weights_functions.R")

##################################################
# 3 - Give values to the following parameters
##################################################

# shape: shape parameter of the Weibull distribution of the baseline hazard function for T
# scale: scale parameter of the Weibull distribution of the baseline hazard function for T
# n.var: number of covariates (normally distributed)
# mean: mean vector of the covariates
# covariance: covariance matrix of the covariates
# beta: vector of coefficients for the cox regression
# n.sample: sample size 
# n.rep: number of partitions
# i: seed

################################################
# 4 - Best model selection
################################################

res<-numeric()
res.alasso<-numeric()
res.alasso.pca<-numeric()
res.alasso.uni<-numeric()
ind<-numeric()
c_index<-c()
c_index1<-c()
c_index2<-c()
c_index3<-c()

dat<-survival.time.weibull(shape=shape,scale=scale,n.sample=n.sample,n.var=n.var, mean=mean, covariance=covariance, beta=beta,i=i)
x_data<-dat$x
y_data<-dat$t
for(i in 1:n.rep){

  index<-sample(1:n.sample,n.sample/2)
  ind<-cbind(ind,index)
  x<-x_data[index,]
  x_test<-x_data[-index,]
  y<-y_data[index,]
  y_test<-y_data[-index,]
  
  #lasso
  cvfit <- cv.glmnet(x, y, family = "cox", alpha=1)
  coeficients<-coef(cvfit,s=cvfit$lambda.min)
  coef<-c(which(coef(cvfit,s=cvfit$lambda.min)!=0))
  data<-as.data.frame(x)
  data$y1<-y[,1]
  data$y2<-y[,2]
  
  data_test<-as.data.frame(x_test)
  data_test$y1<-y_test[,1]
  data_test$y2<-y_test[,2]
  coeff<-c()
  for (k in coef){coeff<-c(coeff,paste("data[,",k,"]"))}
  cox<-coxph(as.formula(paste('Surv(y1,y2) ~', paste(c(coeff), collapse='+'))),data)
  lpnew <- predict(cox, newdata=data_test)
  c_ind1<-GHCI(lpnew)
  c_index<-c(c_index,cind1)
  
  
  #adaptive lasso with ridge weights
  weights_ridge<-Ridge_weights(x,y)
  cvfit.alasso <- cv.glmnet(x, y, family = "cox", alpha=1, penalty.factor=weights_ridge)
  coeficients.alasso<-coef(cvfit.alasso,s=cvfit.alasso$lambda.min)
  coef.alasso<-c(which(coef(cvfit.alasso,s=cvfit.alasso$lambda.min)!=0))
  coeff.alasso<-c()
  for (k in coef.alasso){coeff.alasso<-c(coeff.alasso,paste("data[,",k,"]"))}
  cox.alasso<-coxph(as.formula(paste('Surv(y1,y2) ~', paste(c(coeff.alasso), collapse='+'))),data)
  lpnew.alasso <- predict(cox.alasso, newdata=data_test)
  cind2<-GHCI(lpnew.alasso)
  c_index1<-c(c_index1,cind2)

  
  #adaptive lasso with PCA weights
  weights_PCA<-PCA_weights(x,y)
  cvfit.alasso.pca <- cv.glmnet(x, y, family = "cox", alpha=1, penalty.factor=weights_PCA)
  coeficients.pca<-coef(cvfit.alasso.pca,s=cvfit.alasso.pca$lambda.min)
  coef.pca<-c(which(coef(cvfit.alasso.pca,s=cvfit.alasso.pca$lambda.min)!=0))
  coeff.pca<-c()
  for (k in coef.pca){coeff.pca<-c(coeff.pca,paste("data[,",k,"]"))}
  cox.pca<-coxph(as.formula(paste('Surv(y1,y2) ~', paste(c(coeff.pca), collapse='+'))),data)
  lpnew.pca <- predict(cox.pca, newdata=data_test)
  cind3<-GHCI(lpnew.pca)
  c_index2<-c(c_index2,cind3)

  
  #adaptive lasso with univariate weights
  weights_uni<-uni_weights(x,y)
  cvfit.alasso.uni <- cv.glmnet(x, y, family = "cox", alpha=1, penalty.factor=weights_uni)
  coeficients.uni<-coef(cvfit.alasso.uni,s=cvfit.alasso.uni$lambda.min)
  coef.uni<-c(which(coef(cvfit.alasso.uni,s=cvfit.alasso.uni$lambda.min)!=0))
  coeff.uni<-c()
  for (k in coef.uni){coeff.uni<-c(coeff.uni,paste("data[,",k,"]"))}
  cox.uni<-coxph(as.formula(paste('Surv(y1,y2) ~', paste(c(coeff.uni), collapse='+'))),data)
  lpnew.uni <- predict(cox.uni, newdata=data_test)
  cind4<-GHCI(lpnew.uni)
  c_index3<-c(c_index3,cind4)
  
  
  res<-rbind(res,coeficients[,1])
  res.alasso<-rbind(res.alasso,coeficients.alasso[,1])
  res.alasso.pca<-rbind(res.alasso.pca,coeficients.alasso.pca[,1])
  res.alasso.uni<-rbind(res.alasso.uni,coeficients.alasso.uni[,1])
}

# Best model selection for Lasso
imp<-c()
for (j in 1:ncol(x)){
  importance<-(t(abs(res[,j]))%*%c_index)
  imp<-c(imp,importance)
}
imp_index<-imp/max(imp)
sorted_imp<-sort(imp_index, decreasing=TRUE,index.return=TRUE)$x
sorted_imp_ind<-sort(imp_index, decreasing=TRUE,index.return=TRUE)$ix
k<-ceiling(sqrt(n.sample/2))
power<-c()
for (j in 1:n.rep){
  power[j]<-sum(t(abs(res[j,sorted_imp_ind[k:n.var]]))%*%sorted_imp[k:n.var])/(sum(sorted_imp[1:k])*sum(abs(res[j,])))
}
best_model<-res[which(power+c_index==max(power+c_index0)),res[which(power+c_index==max(power+c_index)),]!=0]
best_model_covariates<-names(best_model)
