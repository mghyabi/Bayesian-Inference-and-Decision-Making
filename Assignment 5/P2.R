if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

require(rmutil)

alpha0=1
beta0=3

priorMean=alpha0/(alpha0+beta0)
priorSD=sqrt((alpha0*beta0)/((alpha0+beta0)^2*(alpha0+beta0+1)))

priorLower = qbeta(0.025, shape1=alpha0, shape2=beta0)  # lower bound on 95% interval
priorUpper = qbeta(0.975, shape1=alpha0, shape2=beta0)  # lower bound on 95% interval

pi=seq(length=100,from=0.01,to=1)
priorDens=dbeta(pi,shape1=alpha0,shape2=beta0)
plot(pi,priorDens,type="l",col="blue",
     main="Prior Distribution",
     xlab=expression(pi),ylab="Probability Density",
     xlim=c(0,max(pi)),ylim=c(0,max(priorDens)))

#----

alpha1=alpha0+19
beta1=beta0+47-19

postMean=alpha1/(alpha1+beta1)
postSD=sqrt((alpha1*beta1)/((alpha1+beta1)^2*(alpha1+beta1+1)))

postLower = qbeta(0.025, shape1=alpha1, shape2=beta1)  # lower bound on 95% interval
postUpper = qbeta(0.975, shape1=alpha1, shape2=beta1)  # lower bound on 95% interval

postDens=dbeta(pi,shape1=alpha1,shape2=beta1)
normLikC=dbeta(pi,shape1=1+19,shape2=1+47-19)
plot(pi,priorDens,type="l",col="blue",
     main="Triplot",
     xlab=expression(pi),ylab="Probability Density",
     xlim=c(0,max(pi)),ylim=c(0,6.5))
lines(pi,normLikC,col="green")
lines(pi,postDens,col="red")
legend(0.8,5.0,c("Prior","Norm Lik","Posterior"),col=c("blue","green","red"),lty=c(1,1,1))

#------
x=0:50

predDist=dbetabinom(x, size = 50, alpha1/(alpha1+beta1), alpha1+beta1)
test=dbinom(x, size=50, prob=19/47)

barplot(predDist,main="Predictive Distribution", 
        xlab="Number of people who choose B and C", 
        ylab="Probability", col="lightblue", 
        border="darkblue",names.arg =x,
        ylim=c(0,0.1))
Data <- rbind(predDist,test)
barplot(Data,main="Predictive Dist. compared with binomial Dist. with sample frequency as point estimate probability", 
        xlab="Number of Observations in a 15-second Block", 
        ylab="EProbability", col=c("lightblue","pink"), 
        border=c("darkblue","red"),names.arg =x, beside=TRUE,
        legend=c("Predictive Distribution","point estimate of binomial"),ylim=c(0,0.15))
test2=dbinom(x, size=50, prob=postMean)
Data <- rbind(predDist,test2)
barplot(Data,main="Predictive Dist. compared with binomial Dist with posterior mean as point estimate probability", 
        xlab="Number of Observations in a 15-second Block", 
        ylab="EProbability", col=c("lightblue","pink"), 
        border=c("darkblue","red"),names.arg =x, beside=TRUE,
        legend=c("Predictive Distribution","point estimate of binomial"),ylim=c(0,0.15))