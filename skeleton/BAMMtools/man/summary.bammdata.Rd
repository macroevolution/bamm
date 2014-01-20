\name{summary.bammdata}
\alias{summary.bammdata}
\title{
Summary of rate shift results from \code{BAMM} analysis.
}
\description{
Summarizes the number of shifts and their locations in the 
maximum shift credibility (MSC) configuration and the posterior 
distribution on the number of shifts.
}
\usage{
S3 method for class \code{bammdata}:
summary(ephy, opt = "product", display = 10)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{opt}{
a character string indicating the optimality criterion to be used in calculation of the maximum shift credibility configuration. May be "product" or "sum"
}
  \item{display}{
an integer for the number of rows of the posterior to display
}
}
\details{
Prints to console the number of posterior samples, the number of 
shifts in the MSC configuration, the locations of those shifts, the optimality
criterion used for calculating the MSC configuration, and the posterior 
distribution on the number of shifts, which is just the fraction of 
samples in the posterior having 0,1,2,...n shifts.
}
\value{
Returns (invisibly) a list with two components:
\item{posterior}{
	the number of posterior samples having 0,1,2,...,n shifts
}
\item{mscshiftnodes}{
	the node(s) where a shift occurs in the maximum shift credibility
	configuration
}
}
\references{
\url{www.bamm-project.org}
}
\author{
Mike Grundler
}
 
\seealso{
\code{\link{maximumShiftCredibility}}
}
\examples{
library(BAMMtools)
data(whales, events.whales)
ephy = getEventData(whales,events.whales,nsamples=100)
summary(ephy)
}

\keyword{ models }
