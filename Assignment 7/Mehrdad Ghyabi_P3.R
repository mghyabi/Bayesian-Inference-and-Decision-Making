#if(!is.null(dev.list())) dev.off()
#cat("\014") 
#rm(list=ls())

Xs <- c(3.74,4.61,4.00,4.67,4.87,5.12,4.52,5.29,5.74,5.48)
Xb <- c(5.44,6.88,5.37,5.44,5.03,6.48,3.89,5.85,6.85,7.16)

n=10

sigma_s = sd(Xs)
sigma_b = sd(Xb)

mu0 = 0
tau0 = Inf

mu1_s=((mu0/tau0^2)+(sum(Xs)/sigma_s^2))/(tau0^-2+(length(Xs)/sigma_s^2))
mu1_b=((mu0/tau0^2)+(sum(Xb)/sigma_b^2))/(tau0^-2+(length(Xb)/sigma_b^2))

tau1_s=(1/tau0^2+length(Xs)/sigma_s^2)^-0.5
tau1_b=(1/tau0^2+length(Xb)/sigma_b^2)^-0.5

SD_s=sqrt(sigma_s^2/n+tau1_s^2)
SD_b=sqrt(sigma_b^2/n+tau1_b^2)

XbarBounds_s=qnorm(c(0.025,0.975), mean = mu1_s, sd = SD_s)
XbarBounds_b=qnorm(c(0.025,0.975), mean = mu1_b, sd = SD_b)

theta=seq(length=1000,from=3,to=8)

marginal_s=dnorm(theta, mean = mu1_s, sd = SD_s)
marginal_b=dnorm(theta, mean = mu1_b, sd = SD_b)

plot(theta,marginal_s,type="l",col="blue",main="Predictive 10 Samples Mean Distributions",
     xlab=expression(bar(X)),ylab="Probability Density",
     xlim=c(3,8),ylim=c(0,2.0),lwd=2)
lines(theta,marginal_b,col="red",lwd=2)
lines(XbarBounds_s,c(0.1,0.1),type="b",col="blue")
lines(XbarBounds_b,c(0.05,0.05),type="b",col="red")
legend(5.5,1.75,c("Surface Sample Mean","Bottom Sample Mean"),col=c("blue","red"),lty=c(1,1))

n=40

sigma_s = sd(Xs)
sigma_b = sd(Xb)

mu0 = 0
tau0 = Inf

mu1_s=((mu0/tau0^2)+(sum(Xs)/sigma_s^2))/(tau0^-2+(length(Xs)/sigma_s^2))
mu1_b=((mu0/tau0^2)+(sum(Xb)/sigma_b^2))/(tau0^-2+(length(Xb)/sigma_b^2))

tau1_s=(1/tau0^2+length(Xs)/sigma_s^2)^-0.5
tau1_b=(1/tau0^2+length(Xb)/sigma_b^2)^-0.5

SD_s2=sqrt(sigma_s^2/n+tau1_s^2)
SD_b2=sqrt(sigma_b^2/n+tau1_b^2)

XbarBounds_s2=qnorm(c(0.025,0.975), mean = mu1_s, sd = SD_s2)
XbarBounds_b2=qnorm(c(0.025,0.975), mean = mu1_b, sd = SD_b2)

theta=seq(length=1000,from=3,to=8)

marginal_s=dnorm(theta, mean = mu1_s, sd = SD_s2)
marginal_b=dnorm(theta, mean = mu1_b, sd = SD_b2)

plot(theta,marginal_s,type="l",col="blue",main="Predictive 40 Samples Mean Distributions",
     xlab=expression(bar(X)),ylab="Probability Density",
     xlim=c(3,8),ylim=c(0,2.0),lwd=2)
lines(theta,marginal_b,col="red",lwd=2)
lines(XbarBounds_s2,c(0.1,0.1),type="b",col="blue")
lines(XbarBounds_b2,c(0.05,0.05),type="b",col="red")
legend(5.5,1.75,c("Surface Sample Mean","Bottom Sample Mean"),col=c("blue","red"),lty=c(1,1))

