# Functions for calculating different proposals of adaptive weights

# x matrix of covariates with n (sample size) rows and p (number of covariates) column
# y matrix with two columns: the survival times and the censoring indicator (0 if censored, 1 if not)

#Function for calculating ridge weights
Ridge_weights<-function(x,y){
  cvfit.ridge <- cv.glmnet(x, y, family = "cox", alpha=0)
  coeficients.ridge<-coef(cvfit.ridge,s=cvfit.ridge$lambda.min)
  weights<-1/abs(coeficients.ridge[,1])
  return(weights)
}

#Function for calculating PCA weights
PCA_weights<-function(x,y){
  PCA<-prcomp(x,scores=TRUE,loadongs=TRUE)
  T<-PCA$x
  P<-PCA$rotation
  std_dev <- PCA$sdev
  pr_var <- std_dev^2
  prop_varex <- pr_var/sum(pr_var)
  sel<-sum(cumsum(prop_varex)<=0.95)+1
  Pdis<-P[,1:sel]
  Tdis<-T[,1:sel]
  cox <- cv.glmnet(Tdis, y, family = "cox", alpha=0)
  beta<-coef(cox,s=cox$lambda.min)
  beta_final<-Pdis%*%beta
  weights<-1/abs(beta_final)
  return(weights)
}


#Function for calculating univariate weights
uni_weights<-function(x,y){
  coef<-c()
  data<-as.data.frame(x)
  data$y1<-y[,1]
  data$y2<-y[,2]
  for (j in 1:4000){ 
    uni<- coxph(Surv(y1,y2)~data[,j],data)
    coef<-c(coef,uni$coefficients[[1]])
  }
  weights<-1/abs(coef)
  return(weights)
}


#Function for calculating random survival forest weights
rsf_weights<-function(x,y){
  data<-as.data.frame(x)
  data$y1<-y[,1]
  data$y2<-y[,2]
  rsf<- rfsrc(Surv(y1,y2)~.,data,importance=TRUE)
  rf.coefs<-rsf$importance
  weights<-1/abs(rf.coefs+(10^-19))
  return(weights)
}

