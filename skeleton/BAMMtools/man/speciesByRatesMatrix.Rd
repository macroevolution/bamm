\name{speciesByRatesMatrix}
\alias{speciesByRatesMatrix}
\title{
	Compute species-specific rate through time trajectories
}
\description{
	Computes the mean of the marginal posterior density of speciation/extinction or phenotypic rates for equally spaced points along the root to tip path for each species
}
\usage{
	speciesByRatesMatrix(ephy, nslices, index, spex)
}
\arguments{
	\item{ephy}{
		an object of class \code{bammdata}
	}
	\item{nslices}{
		an integer number of time slices. This determines the number of equally spaced points in time at which rates are computed for each species
	}
	\item{index}{
		an integer or vector of mode integer indicating which posterior samples to use in the calculation. If \code{NULL} (default) all samples are used
	}
	\item{spex}{
		a character string. "s" (default) calculates speciation rates; "e" calculates extinction rates; "se" calculates diversification rates. Ignored if \code{ephy$type = "trait"}
	}
}

\value{
	A list with two components:
	\item{times}{
		A vector of time points where rates were calculated
	}
	\item{rates}{
		A species X times matrix of rate through time trajectories
	}
}
\references{
	\url{bamm-project.org}
}
\author{
	Mike Grundler
}
\seealso{
	\code{\link{getRateThroughTimeMatrix}}
}
\examples{
	data(whales, events.whales)
	ed <- getEventData(whales,events.whales,burnin=0.25)
	ratemat <- speciesByRatesMatrix(ed, nslices = 100)
	
	dolphins <- extract.clade(whales,140)$tip.label
	plot.new()
	plot.window(xlim=c(0,35),ylim=c(0,0.8))
	for (i in 1:nrow(ratemat$rates)) {
		if (whales$tip.label[i] \%in\% dolphins) {
			lines(ratemat$times, ratemat$rates[i,], lwd=2, col=4)	
		}
		else {
			lines(ratemat$times, ratemat$rates[i,], lwd=2, col=8)
		}
	}
	axis(1,seq(-5,35,5))
	axis(2,seq(-0.2,0.8,0.2),las=1)
	mtext("Time since root",1,line=2.5)
	mtext("Speciation rate",2,line=2.5)
	
}
\keyword{ models }
