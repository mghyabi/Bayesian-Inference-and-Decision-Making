Data <- read.table("http://www.biostat.umn.edu/~lynn/iid/estriol.dat", header=T)

n <- nrow(Data)

x <- Data$estriol
y <- Data$birthwt
xbar <- mean(x)
ybar <- mean(y)

Sxx <- sum((x-xbar)^2)
Syy <- sum((y-ybar)^2)
Sxy <- sum((x-xbar)*(y-ybar))

b=Sxy/Sxx
a=ybar-b*xbar

See <- sum((y-ybar-b*(x-xbar))^2)

rho.shape= (n-2)/2
rho.scale=2/See
rho.bounds=qgamma(c(0.05,0.95), shape = rho.shape, scale = rho.scale)
rho=seq(length=1000,from=0,to=0.2)
plot(rho,dgamma(rho, shape = rho.shape, scale = rho.scale),
     type="l",col="blue",main=bquote("Posterior Distribution of" ~rho),
     xlab=expression(Rho),ylab="Probability Density",
     xlim=c(0,0.2),ylim=c(0,25.0),lwd=2)
lines(rho.bounds,c(0.1,0.1),type="b",col="blue")

eta.center=ybar
eta.spread=1/sqrt(n*(n-2)/See)
eta.df=n-2
eta.bounds=qt(c(0.05,0.95),eta.df)*eta.spread+eta.center
eta=seq(length=1000,from=28,to=36)
plot(eta,dt((eta-eta.center)/eta.spread,eta.df)/eta.spread,
     type="l",col="blue",main=bquote("Posterior Distribution of" ~eta),
     xlab=expression(Eta),ylab="Probability Density",
     xlim=c(min(eta),max(eta)),ylim=c(0,0.65),lwd=2)
lines(eta.bounds,c(0.025,0.025),type="b",col="blue")

beta.center=b
beta.spread=1/sqrt(Sxx*(n-2)/See)
beta.df=n-2
beta.bounds=qt(c(0.05,0.95),beta.df)*beta.spread+beta.center
beta=seq(length=1000,from=0,to=1.2)
plot(beta,dt((beta-beta.center)/beta.spread,beta.df)/beta.spread,
     type="l",col="blue",main=bquote("Posterior Distribution of" ~beta),
     xlab=expression(Beta),ylab="Probability Density",
     xlim=c(min(beta),max(beta)),ylim=c(0,3.0),lwd=2)
lines(beta.bounds,c(0.025,0.025),type="b",col="blue")

#-------

plot(x,y, main = "Estriol vs Birth Weight",
     xlab="Estriol (mg/24h)", ylab="Birth Weight (g/100)")
lines(c(-10,50),a+b*c(-10,50),col="blue",lwd=2)

qqnorm(y-(a+b*x), main = "Residual Q-Q plot")
qqline(y-(a+b*x),col="blue",lwd=2)

#----

xnew=19
pred.center=xnew*b+a
pred.spread= sqrt(((xnew-xbar)^2/Sxx+1/n+1)*(See/(n-2)))
pred.df=n-2
pred.bounds=qt(c(0.05,0.95),pred.df)*pred.spread+pred.center
pred=seq(length=1000,from=0,to=60)
plot(pred,dt((pred-pred.center)/pred.spread,pred.df)/pred.spread,
     type="l",col="blue",main="Birthweight PDF given Estriol=19 mg/24h",
     xlab="Birth Weight (g/100)",ylab="Probability Density",
     xlim=c(min(pred),max(pred)),ylim=c(0,0.12),lwd=2)
lines(pred.bounds,c(0.0025,0.0025),type="b",col="blue")

#----

nsim=1000
xnew=19

rhoMC = etaMC = betaMC = yMC = NULL
for (k in 1:nsim) {
  rhoMC[k]=rgamma(1,shape = (n-2)/2, scale = 2/See)
  etaMC[k]=rnorm(1,mean = ybar, sd=1/sqrt(n*rhoMC[k]))
  betaMC[k]=rnorm(1,mean = b, sd=1/sqrt(Sxx*rhoMC[k]))
  yMC[k]=rnorm(1,mean = etaMC[k]+betaMC[k]*(xnew-xbar), sd=sqrt(1/rhoMC[k]))
}
yMC.bounds=quantile(yMC,c(0.05,0.95))

plot(density(yMC),main="Estimated Birthweight PDF given Estriol=19 mg/24h",
     type="l",col="blue",
     xlab="Birth Weight (g/100)",ylab="Probability Density",
     xlim=c(min(pred),max(pred)),ylim=c(0,0.12),lwd=2)
lines(pred,dt((pred-pred.center)/pred.spread,pred.spread)/pred.spread,
      col="blue",lty=2)
lines(yMC.bounds,c(0.0025,0.0025),type="b",col="blue")
legend(3,0.1,c("MC Sampling","Theoretical"),col=c("blue","blue"),lty=c(1,2))

