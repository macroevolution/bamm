\name{getJenksBreaks}
\alias{getJenksBreaks}
\title{
Jenks natural breaks classification
}
\description{
	Given a vector of numeric values and the number of desired categories, calculate the optimum breakpoints
	using Jenks natural breaks optimization.
}
\usage{
getJenksBreaks(var, k, subset = NULL)
}
\arguments{
 \item{var}{
 	numeric vector
}
\item{k}{
	number of groupings
}
\item{subset}{
	Number of regularly spaced samples to subset from \code{var}. 
	Intended to improve runtime for large datasets.
	If \code{NULL}, all values are used.
}
}

\details{
	\code{getJenksBreaks} is called by \code{\link{assignColorBreaks}}.

	The values in \code{var} are binned into \code{k} categories, according to the Jenks natural 
	breaks classification method. This method is borrowed from the field of cartography, and 
	seeks to minimize the variance within categories, while maximizing the variance between categories. 
	If \code{subset = NULL}, all values of \code{var} are used for the optimization, however this
	can be a slow process with very large datasets. If \code{subset} is set to some number, then 
	\code{subset} regularly spaced values of \code{var} will be sampled. This is slightly less accurate
	than when using the entirety of \code{var} but is unlikely to make much of a difference. If 
	\code{subset} is defined but \code{length(var) < subset}, then \code{subset} has no effect.  


	The Jenks natural breaks method was ported to C from code found in the classInt R package. 

}

\value{
A numeric vector of intervals.
}

\author{
Pascal Title
}

\seealso{
	See \code{\link{assignColorBreaks}} and \code{\link{plot.bammdata}}. 
}

\examples{
#load whales dataset
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)

#for demonstration purposes, we will extract the vector of speciation rates
ed <- dtRates(ed, tau=0.01)
vec <- ed$dtrates$rates[[1]]

#Return breaks for the binning of speciation rates into 64 groups
getJenksBreaks(vec, 64)

}
 
\keyword{ graphics }
