if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

carCounts <- c(3, 5, 7, 3, 3)  # Counts were given as inputs for HW 4

#prior
alpha0=1
beta0=Inf
#poesterior
alpha1=alpha0+sum(carCounts*c(0:4))
beta1=1/(1/beta0 + sum(carCounts))
#predictive
size=alpha1 
prob=1/(1+beta1)

predMass=dnbinom(0:4, size, prob)  #using a negative binomial distribution for predictive distribution

PMF <- c(0.14847768,0.27670841,0.26413075,0.17208518,0.08604259)  #from Assignment 3

Data <- rbind(c(predMass,1-sum(predMass)),c(PMF,1-sum(PMF)))
barplot(Data,main="Predictive Distribution in Comparison with Results from Assignment 3", 
        xlab="Number of Observations in a 15-second Block", 
        ylab="EProbability", col=c("lightblue","pink"), 
        border=c("darkblue","red"),names.arg =c(0:4,"5 or more"), beside=TRUE,
        legend=c("Predictive Distribution","From Assignment 3"),ylim=c(0,0.3))

