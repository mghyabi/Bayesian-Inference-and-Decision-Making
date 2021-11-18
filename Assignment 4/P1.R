if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

Intervals <- c(12,2,6,2,19,5,34,4,1,4,8,7,1,21,6,11,8,28,6,4,5,1,18,9,5,1,21,1,1,5,3,14,5,3,4,5,1,3,16,2)
NumEvent <- ceiling(cumsum(Intervals/15))

Rates <- c()
for (i in min(NumEvent):max(NumEvent)) {
  Rates <- c(Rates,sum(i==NumEvent))
}

alpha0=1; beta0=Inf                  # prior shape and scale
alpha1=alpha0+sum(Rates)             # posterior shape
beta1=1/(1/beta0+length(Rates))      # posterior scale

lambda=seq(length=100,from=0.04,to=4)
priorDens=dgamma(lambda,shape=alpha0,scale=beta0)
postDens=dgamma(lambda,shape=alpha1,scale=beta1)
normLikC=dgamma(lambda,shape=sum(Rates)+1,scale=1/length(Rates))

plot(lambda,postDens,type="l",col="red",main="Posterior Distribution for Traffic Rate",
     xlab="Traffic Rate (number of cars per 15 seconds)",ylab="Probability Density",
     xlim=c(0,4),ylim=c(0,1.5))

Mean1=alpha1*beta1          # posterior mean
SD1=sqrt(alpha1*beta1^2)  # posterior standard deviation
Mode1=(alpha1-1)*beta1      # posterior mode
Median1=qgamma(0.5, shape=alpha1, scale=beta1) # posterior median
Int95=qgamma(c(0.025,0.975), shape=alpha1, scale=beta1) # posterior 95% interval
