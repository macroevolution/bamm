# Jon: Build polygon
buildPolygon <- function(y)	{
	# Jon: got this from Stackoverflow forever ago
	x <- as.numeric(names(y))
	y2 <- rep(y, each=2)
	y2 <- y2[-length(y2)]
	x2 <- rep(x, each=2)[-1]
	x3 <- c(min(x2), x2, max(x2))
	y3 <- c(0, y2, 0)
	return(list(x=x3, y=y3))	
}

# Jon: Analytical prior giving the probability of observing K shifts given a Poisson Rate Prior
prob.k <- function(k, poissonRatePrior=1)	{
	Denom <- (poissonRatePrior + 1)^(k+1)
	Prob <- poissonRatePrior / Denom
	return(Prob)
}

# Jon:Given an mcmc_out output file from BAMM, and the poisson rate prior used, calculate the prior probability of observing 0 to max N shifts
plotPrior <- function(mcmc_out, expectedNumberofShifts, burnin=0.1, plot=TRUE, defs=TRUE, poly=TRUE, pts=FALSE, prCol=rgb(1,0,0,0.5), poCol=rgb(0,0,1,0.5))	{
	Start <- floor(nrow(mcmc_out)*burnin)
	End <- nrow(mcmc_out)
	
	mcmc_out <- mcmc_out[Start:End,]
	
	# Jon: generates a vector from 0 to Limit more than the largest observed number of shifts 
	obsK <- seq(from=0, to=max(mcmc_out[,"N_shifts"]), by=1)
	
	# Jon: calculates the prior probability of 0 to max shifts
	prior <- sapply(obsK, prob.k, poissonRatePrior=1/expectedNumberofShifts)
	
	# Jon: format the output file to be ready to read straight into, e.g., computeBayesFactors()
	out <- data.frame(N_shifts = obsK, prob = prior)
	
	# Jon: optional: plot the prior and posterior against one another. Is pretty.
	if (plot==TRUE)	{
		postD <- tapply(mcmc_out[,"N_shifts"], mcmc_out[,"N_shifts"], length)/nrow(mcmc_out)
		postD <- postD/sum(postD)
		postD <- postD[order(as.numeric(names(postD)))]
		priorD <- out[,2]
		names(priorD) <- out[,1]
		
		# Jon: plotting parameter defaults
		if (defs==TRUE)	{
			par(mar=c(4,4,0,0), mgp=c(2.3,0.5,0), tck=-0.01, las=1, cex.axis=0.8, cex.lab=1.3, bty="n")
		}

		# Jon: make an empty plot...
		plot(1, 1, xlim=c(0,max(obsK)), ylim=c(range(c(postD, priorD))), type="n", xlab="# of shifts", ylab="probability")
		
		# Jon: distros as polygons?
		if (poly==TRUE)	{
			vertices1 <- buildPolygon(priorD)
			vertices2 <- buildPolygon(postD)
			polygon(vertices1$x-0.5, vertices1$y, border=rgb(0,0,0,0), col=prCol)
			polygon(vertices2$x-0.5, vertices2$y, border=rgb(0,0,0,0), col=poCol)
		}
		
		# Jon: distros as points?
		if (pts==TRUE)	{
			points(as.numeric(names(priorD)), priorD, pch=16, col=prCol, type="b")
			points(as.numeric(names(postD)), postD, pch=16, col=poCol, type="b")
		}
	}

	# Jon: return the output file for additional analysis
	return(out)
}
