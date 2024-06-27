# Function for simulating survival times with censoring

# shape: shape parameter for the weibull distribution of the censoring variable C (same value as for the distribution of T)
# theta: scale parameter for the weibull distribution of the censoring variable C. It can be obtained using the function "simulate_scale_censring_proportion.R".
# n.sample: sample size
# T: vector containing the survival times. It can be computed using the function "generate_uncensored_data.R".

censored_data<-function(shape,theta,n.sample,T){
  C<-rweibull(n.sample,shape,theta)
  c<-c(rep(1,n.sample))
  time<-c()
  for (j in 1:n.sample){
    time[j]<-min(t[j,i],C[j])
    if (t[j,i]>C[j]){
      c[j]<-0
    }
  }
  y<-cbind(time,c)
  colnames(y)<-c("time","status")
  y
}
