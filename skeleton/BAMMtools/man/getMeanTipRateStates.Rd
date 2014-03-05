\name{getMeanTipRateStates}
\alias{getMeanTipRateStates}
\title{
Compute tip-specific macroevolutionary rates from \code{bammdata} object
}
\description{
Return speciation, extinction, net diversification, or Brownian motion trait rates for all tips in the phylogeny from \code{BAMM} output.
}
\usage{
getMeanTipRateStates(ephy, ratetype, use.names = TRUE)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{ratetype}{
the type of rate to return: either \code{speciation}, \code{extinction}, \code{netdiv}, or \code{trait}.
}
  \item{use.names}{
logical: should the returned vector be named with taxon names, defaults to \code{TRUE}.
}
}

\value{
Returns a vector of tip rates, named if \code{use.names = TRUE}.
}
\author{
Pascal Title
}
\seealso{
Requires an object of class \code{bammdata} as obtained with \code{\link{getEventData}}. 
}
\examples{
data(whales, events.whales)
ephy <- getEventData(whales, events.whales, burnin=0.25)

#return a vector of net diversification rates at the tips of the phylogeny.
meanrates <- getMeanTipRateStates(ephy, ratetype='netdiv')

# Look at a histogram of species-specific rates:
hist(meanrates, breaks=25)

}
 
\keyword{ models }
