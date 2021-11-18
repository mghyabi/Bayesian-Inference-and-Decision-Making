if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

Xs <- c(3.74,4.61,4.00,4.67,4.87,5.12,4.52,5.29,5.74,5.48)
Xb <- c(5.44,6.88,5.37,5.44,5.03,6.48,3.89,5.85,6.85,7.16)

Xbars=mean(Xs)
Xbarb=mean(Xb)

Ys=sum((Xs-Xbars)^2)
Yb=sum((Xb-Xbarb)^2)

mu0=0
kappa0=0
alpha0=-0.5
beta0=Inf

mu1s=(kappa0*mu0+length(Xs)*Xbars)/(kappa0+length(Xs))
mu1b=(kappa0*mu0+length(Xb)*Xbarb)/(kappa0+length(Xb))
kappa1s=kappa0+length(Xs)
kappa1b=kappa0+length(Xb)
alpha1s=alpha0+length(Xs)/2
alpha1b=alpha0+length(Xb)/2
beta1s=1/(1/beta0+0.5*Ys+0.5*length(Xs)*kappa0*(Xbars-mu0)^2/(length(Xs)+kappa0))
beta1b=1/(1/beta0+0.5*Yb+0.5*length(Xb)*kappa0*(Xbarb-mu0)^2/(length(Xb)+kappa0))

RHObounds=qgamma(c(0.05,0.95), shape = alpha1s, scale = beta1s)
RHOboundb=qgamma(c(0.05,0.95), shape = alpha1b, scale = beta1b)

spreads=(kappa1s*alpha1s*beta1s)^-0.5
spreadb=(kappa1b*alpha1b*beta1b)^-0.5

THETAbounds=c(mu1s-qt(0.95,2*alpha1s)*spreads,mu1s+qt(0.95,2*alpha1s)*spreads)
THETAboundb=c(mu1b-qt(0.95,2*alpha1b)*spreadb,mu1b+qt(0.95,2*alpha1b)*spreadb)

rho=seq(length=500,from=0,to=5)
rhoPostDensS=dgamma(rho,shape=alpha1s,scale=beta1s)  
rhoPostDensB=dgamma(rho,shape=alpha1b,scale=beta1b)  

plot(rho,rhoPostDensS,type="l",col="blue",main="Posterior Precision Distributions",
     xlab=expression(Rho),ylab="Probability Density",
     xlim=c(0,5),ylim=c(0,1.0),lwd=2)
lines(rho,rhoPostDensB,col="red",lwd=2)
legend(3.0,1.0,c("Surface Precision","Bottom Precision"),col=c("blue","red"),lty=c(1,1))

theta=seq(length=500,from=3,to=8)
TS=(theta-mu1s)/spreads
TB=(theta-mu1b)/spreadb
thetaPostDensS=dt(TS,2*alpha1s)/spreads  
thetaPostDensB=dt(TB,2*alpha1b)/spreadb  

plot(theta,thetaPostDensS,type="l",col="blue",main="Posterior Mean Distributions",
     xlab=expression(Theta),ylab="Probability Density",
     xlim=c(3,8),ylim=c(0,2.0),lwd=2)
lines(theta,thetaPostDensB,col="red",lwd=2)
legend(6.5,1.5,c("Surface Mean","Bottom Mean"),col=c("blue","red"),lty=c(1,1))


#------
#Problem2

numSim <- 10000
rhoDirects <- rgamma(numSim,shape=alpha1s,scale=beta1s) 
rhoDirectb <- rgamma(numSim,shape=alpha1b,scale=beta1b) 
thetaDirects <- rnorm(numSim,mean=mu1s,sd=(kappa1s*rhoDirects)^-0.5)
thetaDirectb <- rnorm(numSim,mean=mu1b,sd=(kappa1b*rhoDirectb)^-0.5)

CDFthetas=cumsum(unlist(density(thetaDirects)[2]))
CDFthetab=cumsum(unlist(density(thetaDirectb)[2]))

THETAbounds2=quantile(thetaDirects,c(0.05,0.95))
THETAboundb2=quantile(thetaDirectb,c(0.05,0.95))
RHOAbounds2=quantile(rhoDirects,c(0.05,0.95))
RHOAboundb2=quantile(rhoDirectb,c(0.05,0.95))

plot(density(thetaDirects),col="blue",lty=1,lwd=1,
     main = "Monte Carlo Estimation Compared to t Distribution",
     xlab=expression(Theta),ylab="Probability Density",
     xlim=c(3,8),ylim=c(0,2.0))
lines(theta,thetaPostDensS,col="blue",lty=2)
lines(density(thetaDirectb),col="red",lty=1,lwd=1,main="",xlab="Theta")
lines(theta,thetaPostDensB,col="red",lty=2)
legend(6.50,1.75,c("Surface Monte Carlo","Surface Theoretical t","Bottom Monte Carlo","Bottom Theoretical t"),
       col=c("blue","blue","red","red"),lty=c(1,2,1,2))

plot(density(rhoDirects),col="blue",lty=1,lwd=1,
     main = "Monte Carlo Estimation Compared to Gamma Distribution",
     xlab=expression(Rho),ylab="Probability Density",
     xlim=c(0,5),ylim=c(0,1.0))
lines(rho,rhoPostDensS,col="blue",lty=2)
lines(density(rhoDirectb),col="red",lty=1,lwd=1,main="",xlab="Theta")
lines(rho,rhoPostDensB,col="red",lty=2)
legend(3,0.9,c("Surface Monte Carlo","Surface Theoretical gamma","Bottom Monte Carlo","Bottom Theoretical gamma"),
       col=c("blue","blue","red","red"),lty=c(1,2,1,2))
#------
#Problem3

sigmaDirects <- 1/sqrt(rhoDirects)
sigmaDirectb <- 1/sqrt(rhoDirectb)

thetaDiff=thetaDirectb-thetaDirects
sigmaDiff=sigmaDirectb-sigmaDirects

ProbThetabHigher=sum(thetaDiff>0)/length(thetaDiff)
ProbSigmabHigher=sum(sigmaDiff>0)/length(sigmaDiff)

#------
#Problem3

qqnorm(Xb,main = "Q-Q plot for Bottom Readings")
qqline(Xb,lty=2,col="blue")
qqnorm(Xs,main = "Q-Q plot for Surface Readings")
qqline(Xs,lty=2,col="blue")