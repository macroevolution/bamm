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
\details{
	
}
\value{
	A nspecies x nslices matrix of rate through time trajectories.
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
	ratemat<-speciesByRatesMatrix(ed, nslices = 50)
	
	pca <- princomp(ratemat)
	plot(pca$scores[,1], pca$scores[,2])
}
\keyword{ models }
