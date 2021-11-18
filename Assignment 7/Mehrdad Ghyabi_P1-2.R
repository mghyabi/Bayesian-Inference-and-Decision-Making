#Problem 1

#if(!is.null(dev.list())) dev.off()
#cat("\014") 
#rm(list=ls())

Xs <- c(3.74,4.61,4.00,4.67,4.87,5.12,4.52,5.29,5.74,5.48)
Xb <- c(5.44,6.88,5.37,5.44,5.03,6.48,3.89,5.85,6.85,7.16)

n=10

mu0 = 0
k0 = 0
alpha0 = -0.5
beta0 = Inf

Xbar_s=mean(Xs)
Xbar_b=mean(Xb)

Ys=sum((Xs-Xbar_s)^2)
Yb=sum((Xb-Xbar_b)^2)

mu1_s=(k0*mu0+length(Xs)*Xbar_s)/(k0+length(Xs))
mu1_b=(k0*mu0+length(Xb)*Xbar_b)/(k0+length(Xb))

k1_s=k0+length(Xs)
k1_b=k0+length(Xb)

alpha1_s=alpha0+length(Xs)/2
alpha1_b=alpha0+length(Xb)/2

beta1_s=1/(1/beta0+0.5*Ys+0.5*length(Xs)*k0*(Xbar_s-mu0)^2/(length(Xs)+k0))
beta1_b=1/(1/beta0+0.5*Yb+0.5*length(Xb)*k0*(Xbar_b-mu0)^2/(length(Xb)+k0))

spread_sP=(k1_s*n*alpha1_s*beta1_s/(k1_s+n))^-0.5
spread_bP=(k1_b*n*alpha1_b*beta1_b/(k1_b+n))^-0.5

XbarBounds_s=c(mu1_s-qt(0.975,2*alpha1_s)*spread_sP,mu1_s+qt(0.975,2*alpha1_s)*spread_sP)
XbarBounds_b=c(mu1_b-qt(0.975,2*alpha1_b)*spread_bP,mu1_b+qt(0.975,2*alpha1_b)*spread_bP)

theta=seq(length=1000,from=3,to=8)
TS=(theta-mu1_s)/spread_sP
TB=(theta-mu1_b)/spread_bP
marginal_s=dt(TS,2*alpha1_s)/spread_sP  
marginal_b=dt(TB,2*alpha1_b)/spread_bP

plot(theta,marginal_s,type="l",col="blue",main="Predictive 40 Samples Mean Distributions",
     xlab=expression(bar(X)),ylab="Probability Density",
     xlim=c(3,8),ylim=c(0,1.75),lwd=2)
lines(theta,marginal_b,col="red",lwd=2)
lines(XbarBounds_s,c(0.1,0.1),type="b",col="blue")
lines(XbarBounds_b,c(0.05,0.05),type="b",col="red")
legend(5.5,1.5,c("Surface Sample Mean","Bottom Sample Mean"),col=c("blue","red"),lty=c(1,1))

#-----
#Problem 2
numSim=10000

XbarDirect_s=rt(numSim,2*alpha1_s)*spread_sP+mu1_s  
XbarDirect_b=rt(numSim,2*alpha1_b)*spread_bP+mu1_b

XbarBounds_s2=quantile(XbarDirect_s,c(0.025,0.975))
XbarBounds_b2=quantile(XbarDirect_b,c(0.025,0.975))

diffBounds=quantile(XbarDirect_b-XbarDirect_s,c(0.025,0.975))

plot(density(XbarDirect_s),col="blue",lty=1,lwd=1,
     main = "Monte Carlo 40 Samples Mean Distributions",
     xlab=expression(bar(X)),ylab="Probability Density",
     xlim=c(3,8),ylim=c(0,2.0))
lines(theta,marginal_s,col="blue",lty=2)
lines(density(XbarDirect_b),col="red",lty=1,lwd=1,main="",xlab="Theta")
lines(theta,marginal_b,col="red",lty=2)
lines(XbarBounds_s2,c(0.1,0.1),type="b",col="blue")
lines(XbarBounds_b2,c(0.05,0.05),type="b",col="red")
legend(5.50,1.75,c("Surface MC Sampling","Surface Theoretical","Bottom MC Sampling","Bottom Theoretical"),
       col=c("blue","blue","red","red"),lty=c(1,2,1,2))

plot(density(XbarDirect_b-XbarDirect_s),col="blue",lty=1,lwd=2,
     main = "Density Function for The Difference in Means for 40 Samples",
     xlab=expression(Theta[b]-Theta[s]),ylab="Probability Density",
     xlim=c(-2,4),ylim=c(0,1.0))
lines(diffBounds,c(0.05,0.05),type="b",col="blue")