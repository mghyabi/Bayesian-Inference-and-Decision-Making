

ReacTimes <- read.table (paste(getwd(),"/Data.txt", sep=""), header=F)
LogTimes <- log(ReacTimes)
LogTimes <- transform(LogTimes, Mean=apply(LogTimes[1:30],1, mean, na.rm = TRUE))
LogTimes <- transform(LogTimes, SD=apply(LogTimes[1:30],1, sd, na.rm = TRUE))
Sigma=mean(LogTimes$SD)
Mean=mean(LogTimes$Mean)

Ind=rep(seq(from=1,to=11),30)
ts=unlist(LogTimes[1:30], use.names = FALSE)
nI=11
n=11*30

IndMeans=LogTimes$Mean
IndSD=LogTimes$SD
IndSize=30

library(rjags)
library(R2jags)

TimeData=list("nI","n","Ind","ts")
TimeParams=c("mu","lambda","theta")
TimeInits <- function(){
  list("mu"=c(5.25),"lambda"=c(25),"theta"=array(5.25,nI))
}

numSim=5000

TimeFit <- jags(data=TimeData, inits=TimeInits, TimeParams, n.chains=1,
                n.iter=numSim, n.burnin=0, model.file="JAGSmodel.jags",
                n.thin=1)
lambdaMC=TimeFit$BUGSoutput$sims.list$lambda
tauMC=sqrt(1/lambdaMC)
muMC=TimeFit$BUGSoutput$sims.list$mu
thetaMC=TimeFit$BUGSoutput$sims.list$theta

muBounds=quantile(muMC,c(0.025,0.975))
tauBounds=quantile(tauMC,c(0.025,0.975))

plot(density(muMC),col="blue",lty=1,lwd=2,
     main = bquote("Estimated Distribution for"  ~mu),
     xlab=expression(mu),ylab="Probability Density",
     xlim=c(5.4,6.0),ylim=c(0,12.0))
lines(muBounds,c(0.2,0.2),type="b",col="blue")
plot(density(tauMC),col="blue",lty=1,lwd=2,
     main = bquote("Estimated Distribution for"  ~tau),
     xlab=expression(tau),ylab="Probability Density",
     xlim=c(0,0.4),ylim=c(0,17.0))
lines(tauBounds,c(0.2,0.2),type="b",col="blue")

q025 <- NULL
q975 <- NULL
for(i in 1:nI) {
  q025[i]<-quantile(thetaMC[,i],0.025)
  q975[i]<-quantile(thetaMC[,i],0.975)
}

require(coda)
traceplot(as.mcmc(muMC), main=bquote("Traceplot for Gibbs Sampling Estimate of"  ~mu))
acf(muMC,main=bquote("ACF for Gibbs Sampling Estimate of"  ~mu))
traceplot(as.mcmc(tauMC), main=bquote("Traceplot for Gibbs Sampling Estimate of"  ~tau))
acf(tauMC,main=bquote("ACF for Gibbs Sampling Estimate of"  ~tau))

#Problem4------
OBS=unlist(LogTimes[,1:30])
qqnorm(OBS,main = "Q-Q plot for all the readings")
qqline(OBS,lty=2,col="blue")


library(MCMCvis)
library(histogram)
TimeFit.mcmc <- as.mcmc(TimeFit)        # change to mcmc object
MCMCsummary(TimeFit.mcmc, round=2)       # show summary of JAGS object

paramMeans=MCMCpstr(TimeFit,          # find parameter means
                    func = mean,
                    type = 'summary')

# Extract effective sample sizes for all variables
effSiz = MCMCpstr(TimeFit,          # find parameter means
                  func = effectiveSize,
                  type = 'summary')
histogram(effSiz$theta,xlab="Effective Sample Size for School Means")

# Shrinkage plots
par(mfrow=c(1,1))
plot(LogTimes$Mean,paramMeans$theta,
     ylab=expression(hat(theta)),xlab=expression(bar(y)))
lines(LogTimes$Mean,LogTimes$Mean)

plot(rep(30,11),(LogTimes$Mean-paramMeans$theta),
     ylab=expression(bar(y)-hat(theta)),xlab="Samples from each individual")
lines(c(0,50), c(0,0))

# Use MCMCvis package to do caterpillar plot for means 
MCMCplot(TimeFit,              
         params = 'theta',
         ref_ovl = TRUE,
         ref=Mean,
         horiz = FALSE
)

MCMCtrace(TimeFit)
