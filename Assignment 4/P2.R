if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

Pars <- get.gamma.par(p = c(0.05, 0.5, 0.95), q = c(1.5,3.75,7), show.output = FALSE, plot = FALSE)

Check <- qgamma(c(0.1,0.50,0.9), shape=Pars[1], scale=1/Pars[2]) 