# Function for calculating the value of the scale parameter for the Weibull distribution of the censoring variable C for obtaining a desired censoring proportion p.

# shape: shape parameter of the Weibull distribution of the baseline hazard function for T
# scale: scale parameter of the Weibull distribution of the baseline hazard function for T
# covariance: covariance matrix of the covariates
# beta: vector of coefficients for the cox regression
# p: censoring proportion

censor.theta.weibull<-function(shape,scale,covariance,beta,p){
  # Add the beta_0 term to the vector of coefficients
  beta1<-c(log(scale^(-shape)),beta)
  # Approximate the function g (censoring rate) using the strong law of large numbers
  g<-function(theta,shape,covariance,beta1,n=10000000){
    u<-rlnorm(n,-beta1[1]/shape,sqrt(t(beta1[-1]/shape)%*%covariance%*%(beta1[-1]/shape)))
    mean(1/(1+((theta/u)^shape)))
  }
  # Find the value of theta that gives gamma(theta)=g(theta)-p=0
  theta<-uniroot(function(u,shape,covariance,beta1,p){g(u,shape,covariance,beta1)-p},shape=shape,covariance=covariance,beta1=beta1,p=p, c(0.0000001,500), tol=0.00000001,extendInt = "yes")$root
  theta
}
