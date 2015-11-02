library(BAMMtools)

prob.k <- function(k, poissonRatePrior=1)	{
	Denom <- (poissonRatePrior + 1)^(k+1)
	Prob <- poissonRatePrior / Denom
	return(Prob)
}

png(file = 'rateshifts_prior.png', width=450, height=450);

obsK <- seq(from=0, to=19, by=1)
expectedNumberofShifts <- 1

priorD <- sapply(obsK, prob.k, poissonRatePrior=1/expectedNumberofShifts)

par(mar=c(5,5,1,1))
plot(1, 1, xlim=c(0,max(obsK)), ylim=c(0,0.5), type='n', xlab="Number of shifts", ylab="Probability")
points(obsK, priorD, pch=21, bg='red', type='b')

dev.off();

##### 

png(file = 'rateshifts_prior_0.1.png', width=450, height=450);


obsK <- seq(from=0, to=19, by=1)
expectedNumberofShifts <- 10

priorD <- sapply(obsK, prob.k, poissonRatePrior=1/expectedNumberofShifts)

par(mar=c(5,5,1,1))
plot(1, 1, xlim=c(0,max(obsK)), ylim=c(0,0.5), type='n', xlab="Number of shifts", ylab="Probability")
points(obsK, priorD, pch=21, bg='red', type='b')

dev.off();



