if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())


Intervals <- c(12,2,6,2,19,5,34,4,1,4,8,7,1,21,6,11,8,28,6,4,5,1,18,9,5,1,21,1,1,5,3,14,5,3,4,5,1,3,16,2)

ObservedCounts <- c()
for (i in min(Intervals):max(Intervals)) {
  ObservedCounts <- c(ObservedCounts,sum(i==Intervals))
}

ExpectedCounts=dexp(min(Intervals):max(Intervals),mean(1/Intervals))*length(Intervals)  

IntervalsExp <- rbind(ObservedCounts,ExpectedCounts)
barplot(IntervalsExp,main="Distribution of Intervals between Observations", 
        xlab="Time between Two Observations (s)", 
        ylab="Empirical Count", col=c("lightblue","pink"), 
        border=c("darkblue","red"),names=min(Intervals):max(Intervals), beside=TRUE,
        legend=c("Observed","Expected"))

exponential.quantiles = qexp(ppoints(length(Intervals)))
qqplot(exponential.quantiles, Intervals,main="Exponential Q-Q Plot of Observation Intervals",
       xlab = "Theoretical Exponential Quantiles", ylab = "Empirical Quantiles")
lines(exponential.quantiles,exponential.quantiles*mean(Intervals),lty=2,col="blue")

#-----------------------------------

NumEvent <- ceiling(cumsum(Intervals/15))

Rates <- c()
for (i in min(NumEvent):max(NumEvent)) {
  Rates <- c(Rates,sum(i==NumEvent))
}

ObservedCounts <- c()
for (i in min(Rates):max(Rates)) {
  ObservedCounts <- c(ObservedCounts,sum(i==Rates))
}

ExpectedCounts=dpois(min(Rates):max(Rates),mean(Rates))*length(Rates)  

RatesExp <- rbind(ObservedCounts,ExpectedCounts) 
barplot(RatesExp,main="Distribution of Observations per 15-second Blocks", 
        xlab="Observations in 15-second Block", 
        ylab="Empirical Count", col=c("lightblue","pink"), 
        border=c("darkblue","red"),names=min(Rates):max(Rates), beside=TRUE,
        legend=c("Observed","Expected"))

poisson.quantiles = qpois(ppoints(length(Intervals)),mean(Rates))
qqplot(poisson.quantiles, Rates,main="Poisson Q-Q Plot of Rate of Observations",
       xlab = "Theoretical Exponential Quantiles", ylab = "Empirical Quantiles")
lines(exponential.quantiles,exponential.quantiles,lty=2,col="blue")

#----------------------------------------------

lambda <- seq(length=20, from=0.2, to=4)

priorDist <- replicate(20,1/20)

barplot(priorDist,main="Prior Distribution",
        xlab="Number of Observations in a 15-second Block", ylab="Probability",names.arg=lambda,
        border="darkblue", col="lightblue",ylim=c(0,0.25))

lik <- array(1,length(lambda))
for (i in 1:10) {
  lik <- lik*dpois(Rates[i],lambda)
}

postDist <- priorDist*lik
postDist <- postDist/sum(postDist)

barplot(postDist,main="Poterior Distribution after observing ten 15-second Blocks",
        xlab="Number of Observations in a 15-second Block", ylab="Probability",names.arg=lambda,
        border="darkblue", col="lightblue",ylim=c(0,0.25))

#Features of the posterior distribution
postMean1 = sum(lambda*postDist)
postVar1 = sum((lambda-postMean1)^2*postDist)
postSD1 = sqrt(postVar1)
postMedian1=lambda[sum(cumsum(postDist)<=0.5)]
percent951=lambda[sum(cumsum(postDist)<=0.95)]

#----------------------------------------------

lambda <- seq(length=20, from=0.2, to=4)

priorDist <- postDist


lik <- array(1,length(lambda))
for (i in 11:21) {
  lik <- lik*dpois(Rates[i],lambda)
}

postDist <- priorDist*lik
postDist <- postDist/sum(postDist)

barplot(postDist,main="Poterior Distribution After Observing All the 15-second Blocks",
        xlab="Number of Observations in a 15-second Block", ylab="Probability",names.arg=lambda,
        border="darkblue", col="lightblue",ylim=c(0,0.25))

#Features of the posterior distribution
postMean2 = sum(lambda*postDist)
postVar2 = sum((lambda-postMean2)^2*postDist)
postSD2 = sqrt(postVar2)
postMedian2=lambda[sum(cumsum(postDist)<=0.5)]
percent952=lambda[sum(cumsum(postDist)<=0.95)]

#---------------------------------

X <- c(0,1,2,3,4)

PMF <- c()
for (i in 1:length(X)) {
  PMF <- c(PMF,sum(dpois(X[i],lambda)*postDist))
}

P4=1-sum(PMF)

barplot(c(PMF,P4),main="Probability Mass for the Number of Observations in the Next 15-second Block",
        xlab="Number of Observations in a 15-second Block", ylab="Probability",names.arg=c(X,">4"),
        border="darkblue", col="lightblue",ylim=c(0,0.3))
