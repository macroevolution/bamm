\name{subsetEventData}
\alias{subsetEventData}

\title{
Subset a \code{bammdata} object
}
\description{
Subsets a \code{bammdata} object. Returns a \code{bammdata} object after extracting a specified set of samples from the posterior.
}
\usage{
subsetEventData(ephy, index)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{index}{
A vector of integers corresponding to samples to be extracted from the posterior distribution of shift configurations included in the \code{bammdata} object. 
}
}
\details{
This will result in an error if you attempt to access samples that do not exist in the \code{ephy} data object. For example, if your \code{bammdata} object includes 100 samples from a posterior distribution sampled with \code{BAMM}, you can only attempt to subset with index values 1:100. 
 

}
\value{
A \code{bammdata} object, but only containing \code{index} samples from the posterior in the original object. 
}
 
\author{
Dan Rabosky
}
 
\examples{
data(whales, events.whales)
ed <- getEventData(whales, events.whales)
ed2 <- subsetEventData(ed, index=1)
plot.bammdata(ed2) 
addBAMMshifts(ed2, cex=2)

}
\keyword{ manip}
