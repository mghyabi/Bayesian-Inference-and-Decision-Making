
Xs <- c(3.74,4.61,4.00,4.67,4.87,5.12,4.52,5.29,5.74,5.48)
Xb <- c(5.44,6.88,5.37,5.44,5.03,6.48,3.89,5.85,6.85,7.16)

n=c(length(Xs),length(Xb))
xbar=c(mean(Xs),mean(Xb))

mu0=6
tau0=1.5
alpha0=4.5
beta0=0.19

#-------
#Problem 2

numSim=10000

rGprev <- c(1/var(Xs),1/var(Xb))
rhoGibbs_s = thetaGibbs_s = sigmaGibbs_s = NULL
rhoGibbs_b = thetaGibbs_b = sigmaGibbs_b = NULL
alphaG=alpha0+n/2
for (k in 1:numSim) {
  muG = (mu0/tau0^2 + n*rGprev*xbar) / (1/tau0^2 + n*rGprev)
  tauG = 1/sqrt(1/tau0^2 + n*rGprev)
  thetaGibbs_s[k] <- rnorm(1,mean=muG[1],sd=tauG[1])
  thetaGibbs_b[k] <- rnorm(1,mean=muG[2],sd=tauG[2])
  betaG<-c(1/(1/beta0 + 0.5*sum((Xs-thetaGibbs_s[k])^2)),1/(1/beta0 + 0.5*sum((Xb-thetaGibbs_b[k])^2)))
  rhoGibbs_s[k]<-rgamma(1,shape=alphaG[1],scale=betaG[1])
  rhoGibbs_b[k]<-rgamma(1,shape=alphaG[2],scale=betaG[2])
  sigmaGibbs_s[k]<-1/sqrt(rhoGibbs_s[k])                # calculate new value of sigma
  sigmaGibbs_b[k]<-1/sqrt(rhoGibbs_b[k])                # calculate new value of sigma
  rGprev = c(rhoGibbs_s[k],rhoGibbs_b[k])    # previous value of rho
}

Dif=thetaGibbs_b-thetaGibbs_s

sigmaBounds_s=quantile(sigmaGibbs_s,c(0.05,0.95))
sigmaBounds_b=quantile(sigmaGibbs_b,c(0.05,0.95))
thetaBounds_s=quantile(thetaGibbs_s,c(0.05,0.95))
thetaBounds_b=quantile(thetaGibbs_b,c(0.05,0.95))
DifBounds=quantile(Dif,c(0.05,0.95))

plot(density(thetaGibbs_s),col="blue",lty=1,lwd=2,
     main = "Gibbs Samples Mean Distributions",
     xlab=expression(Theta),ylab="Probability Density",
     xlim=c(3.5,7.5),ylim=c(0,2.0))
lines(density(thetaGibbs_b),col="red",lty=1,lwd=2,main="",xlab="Theta")
lines(thetaBounds_s,c(0.1,0.1),type="b",col="blue")
lines(thetaBounds_b,c(0.05,0.05),type="b",col="red")
legend(5.50,1.75,c("Surface Gibbs Sampling","Bottom Gibbs Sampling"),col=c("blue","red"),lty=c(1,1))

plot(density(Dif),col="blue",lty=1,lwd=2,
     main = "Difference in Mean Samples Distributions",
     xlab=expression(Theta[b]-Theta[s]),ylab="Probability Density",
     xlim=c(-1.0,3),ylim=c(0,1.1))
lines(DifBounds,c(0.05,0.05),type="b",col="blue")

plot(density(sigmaGibbs_s),col="blue",lty=1,lwd=2,
     main = "Gibbs Samples SD Distributions",
     xlab=expression(Rho),ylab="Probability Density",
     xlim=c(0,2.25),ylim=c(0,3.0))
lines(density(sigmaGibbs_b),col="red",lty=1,lwd=2,main="",xlab="Theta")
lines(sigmaBounds_s,c(0.1,0.1),type="b",col="blue")
lines(sigmaBounds_b,c(0.05,0.05),type="b",col="red")
legend(1.25,2.5,c("Surface Gibbs Sampling","Bottom Gibbs Sampling"),col=c("blue","red"),lty=c(1,1))

#-----
# Problem3
require(coda) 

traceplot(as.mcmc(Dif), main=bquote("Traceplot for Gibbs Sampling Estimate of"  ~Theta[b]-Theta[s]))
acf(Dif,main=bquote("ACF for Gibbs Sampling Estimate of"  ~Theta[b]-Theta[s]))
effectiveSize(Dif)
