setwd("C:/Users/mghyabi/Desktop/PhD/Courses/George Mason/SYST 664/Homework/HW9")
ReacTimes <- read.table ("Data.txt", header=F)
LogTimes <- log(ReacTimes)
LogTimes <- transform(LogTimes, Mean=apply(LogTimes[1:30],1, mean, na.rm = TRUE))
LogTimes <- transform(LogTimes, SD=apply(LogTimes[1:30],1, sd, na.rm = TRUE))
tau=sd(LogTimes$Mean)
mu=mean(LogTimes$Mean)
sigma=mean(LogTimes$SD)

SD1=1/sqrt(1/tau^2+30/sigma^2)
mean1=(LogTimes$Mean*30/sigma^2+mu/tau^2)/(30/sigma^2+1/tau^2)

q025 <- NULL
q975 <- NULL
for(i in 1:11) {
  q025[i]<-qnorm(0.025, mean = mean1[i],sd=SD1)
  q975[i]<-qnorm(0.975, mean = mean1[i],sd=SD1)
}