\name{subsetEventData}
\alias{subsetEventData}

\title{
Subset a \code{bammdata} object
}
\description{
Subsets a \code{bammdata} or a \code{credibleshiftset} object. Returns a \code{bammdata} object after extracting a specified set of samples from the posterior.
}
\usage{
subsetEventData(ephy, index)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata} or class \code{credibleshiftset}
}
  \item{index}{
A vector of integers corresponding to samples to be extracted from the posterior distribution of shift configurations included in the \code{bammdata} object. 
}
}
\details{
This will result in an error if you attempt to access samples that do not exist in the \code{ephy} data object. For example, if your \code{bammdata} object includes 100 samples from a posterior distribution sampled with \code{BAMM}, you can only attempt to subset with index values 1:100. 

This function can also be used to access ranked samples from the credible set of rate shift configurations. See examples.

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

# Using the function to pull out distinct shift configs
# from the credible set
data(prior.whales)
priordist <- getBranchShiftPriors(whales, prior.whales)

css <- credibleShiftSet(ed, threshold=priordist)

# Pull out the 2nd most probable shift configuration:
secondbest <- subsetEventData(css, index=2)

# Plotting it and adding rate shift locations:
plot.bammdata(secondbest, lwd=2)
addBAMMshifts(secondbest, cex=2)

# and here is the overall best:
best <- subsetEventData(css, index=1)

# Setting index = k pulls out 
#	the k'th most-probable shift config

}
\keyword{ manip}
