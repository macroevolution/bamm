\name{getShiftNodesFromIndex}
\alias{getShiftNodesFromIndex}

\title{
Identify nodes associated with rate shifts from \code{bammdata} object
}

\description{
Find the node numbers associated with rate shifts for a specified sample from the posterior distribution contained in a \code{bammdata} object.
}
\usage{
getShiftNodesFromIndex(ephy, index)
}

\arguments{
  \item{ephy}{
a \code{bammdata} object
}
  \item{index}{
The index value of the posterior sample from which you want to identify shift nodes. This is \emph{not} the same as the actual generation number of the MCMC sample. If your \code{bammdata} object contains 100 samples from the posterior distribution, the value of \code{index} must range from 1 to 100.
}
}

\value{
A vector of nodes (excluding the root) that define branches on which shifts occurred for the specified sample from the posterior. Will return a numeric of length 0 if no non-root shifts occur in the specified sample.
}
\author{
Dan Rabosky
}

\seealso{
\link{addBAMMshifts}, \link{plot.bammdata}, \link{maximumShiftCredibility}
}
\examples{
data(whales, events.whales)
ed <- getEventData(whales, events.whales)

# Get the maximum shift credibility configuration:
msc <- maximumShiftCredibility(ed)

# GEt the nodes at which shifts occurred in the 
# maximum shift credibility configuration:

getShiftNodesFromIndex(ed, index=msc$sampleindex)

}
\keyword{ models }
