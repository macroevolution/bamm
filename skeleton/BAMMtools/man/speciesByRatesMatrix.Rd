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
		an integer number of time slices. This determines at how many points in time rates are computed for each species
	}
	\item{index}{
		an integer or vector of mode integer indicating which posterior samples to use in the calculation. If \code{NULL} all samples are used
	}
	\item{spex}{
		a character string. "s" (default) calculates speciation rates; "e" calculates extinction rates; "se" calculates diversification rates. Ignored if \code{ephy$type = "trait"}
	}
\value{
	A nspecies x nslices matrix
}
\references{
	bamm-project.org
}
\author{
	Mike Grundler
}
\examples{
	data(whales, events.whales)
	ed <- getEventData(whales,events.whales,burnin=0.25)
	ratemat<-speciesByRatesMatrix(ed, nslices = 50)
	
	pca <- princomp(ratemat)
	plot(pca$scores[,1], pca$scores[,2])
}
\keyword{ models }
}