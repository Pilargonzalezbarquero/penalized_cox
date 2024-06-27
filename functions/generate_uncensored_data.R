# Function for generating uncensored data given that the baseline hazard function comes from a Weibull distribution 

# shape: shape parameter of the Weibull distribution
# scale: scale parameter of the Weibull distribution
# n.sample: sample size
# n.var: number of covariates (normally distributed)
# mean: mean vector of the covariates
# covariance: covariance matrix of the covariates
# beta: vector of coefficients for the cox regression
# i: seed

survival.time.weibull<-function(shape,scale,n.sample,n.var,mean,covariance,beta,i){
  set.seed(i)
  x<-rmvnorm(n.sample,mean=mean,sigma=covariance)
  t<-scale*(-log(1-runif(n.sample))/exp(as.vector(x%*%beta)))^(1/shape)
  set.seed(NULL)
  y<-cbind(t,1)
  colnames(y)<-c("time","status")
  list(x=x,y=y)
}
