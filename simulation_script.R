########################################################################################################
## Script that generates synthetic data and fits the cox regression models on the train set
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
# p: censoring proportion (if censoring==TRUE) 
# n.sample: sample size for the train and test sets
# n.rep: number of iterations

################################################
# 4 - Simulate train sets and fit Cox models
################################################

res<-numeric()
res.alasso<-numeric()
res.alasso.pca<-numeric()
res.alasso.uni<-numeric()
train_x<-array(0,c(n.sample,n.var,n.rep))
train_y<-array(0,c(n.sample,2,n.rep)) 

if (censoring==TRUE){
  theta<-censor.theta.weibull(shape,scale,covariance,beta,p)
    }
for(i in 1:n.rep){
  dat<-survival.time.weibull(shape=shape,scale=scale,n.sample=n.sample,n.var=n.var, mean=mean, covariance=covariance, beta=beta,i=i)
  x<-dat$x
  y<-dat$t
  train_x[,,i]<-x
  train_y[,,i]<-y
  
  if (censoring==TRUE){
    y<-censored_data(shape,theta,n.sample,y[,1])
  }
  
  #lasso
  cvfit <- cv.glmnet(x, y, family = "cox", alpha=1)
  coeficients<-coef(cvfit,s=cvfit$lambda.min)
  
  
  #adaptive lasso with ridge
  weights_ridge<-Ridge_weights(x,y)
  cvfit.alasso <- cv.glmnet(x, y, family = "cox", alpha=1, penalty.factor=weights_ridge)
  coeficients.alasso<-coef(cvfit.alasso,s=cvfit.alasso$lambda.min)
  
  #adaptive lasso with PCA
  weights_PCA<-PCA_weights(x,y)
  cvfit.alasso.pca <- cv.glmnet(x, y, family = "cox", alpha=1, penalty.factor=weights_PCA)
  coeficients.alasso.pca<-coef(cvfit.alasso.pca,s=cvfit.alasso.pca$lambda.min)
  
  #adaptive lasso with univariate weights
  weights_uni<-uni_weights(x,y)
  cvfit.alasso.uni <- cv.glmnet(x, y, family = "cox", alpha=1, penalty.factor=weights_uni)
  coeficients.alasso.uni<-coef(cvfit.alasso.uni,s=cvfit.alasso.uni$lambda.min)
  
  
  res<-rbind(res,coeficients[,1])
  res.alasso<-rbind(res.alasso,coeficients.alasso[,1])
  res.alasso.pca<-rbind(res.alasso.pca,coeficients.alasso.pca[,1])
  res.alasso.uni<-rbind(res.alasso.uni,coeficients.alasso.uni[,1])
}

# TPR and FPR calculation  
res<-cbind(res,apply(res,1,function(x) sum(x[beta_ind]!=0)/length(beta_ind)),apply(res,1,function(x) sum(x[-beta_ind]!=0)/length(x[-beta_ind])),apply(res,1,function(x) sum(x[beta_ind]==0)/length(beta_ind)))
res.alasso<-cbind(res.alasso,apply(res.alasso,1,function(x) sum(x[beta_ind]!=0)/length(beta_ind)),apply(res.alasso,1,function(x) sum(x[-beta_ind]!=0)/length(x[-beta_ind])),apply(res.alasso,1,function(x) sum(x[beta_ind]==0)/length(beta_ind)))
res.alasso.pca<-cbind(res.alasso.pca,apply(res.alasso.pca,1,function(x) sum(x[beta_ind]!=0)/length(beta_ind)),apply(res.alasso.pca,1,function(x) sum(x[-beta_ind]!=0)/length(x[-beta_ind])),apply(res.alasso.pca,1,function(x) sum(x[beta_ind]==0)/length(beta_ind)))
res.alasso.uni<-cbind(res.alasso.uni,apply(res.alasso.uni,1,function(x) sum(x[beta_ind]!=0)/length(beta_ind)),apply(res.alasso.uni,1,function(x) sum(x[-beta_ind]!=0)/length(x[-beta_ind])),apply(res.alasso.uni,1,function(x) sum(x[beta_ind]==0)/length(beta_ind)))

 ################################################
# 4 - Simulate test sets
#################################################                                                                                                                                                                                                 

test_x<-array(0,c(n.sample,n.var,n.rep))
test_y<-array(0,c(n.sample,2,n.rep))  
                                                                                                                                                                                                  
for(i in 1:n.rep){
  dat<-survival.time.weibull(shape=shape,scale=scale,n.sample=n,n.var=n.var, mean=mean, covariance=covariance, beta=beta,i=i,theta=theta)
  x<-dat$x
  y<-dat$t
  if (censoring==TRUE){
    y<-censored_data(shape,theta,n.sample,y[,1])
  }
  test_x[,,i]<-x
  test_y[,,i]<-y
}

