\name{summary.bammdata}
\alias{summary.bammdata}
\title{
Summary of rate shift results from \code{BAMM} analysis.
}
\description{
Summarizes the posterior distribution on the number of shifts.
}
\usage{
S3 method for class \code{bammdata}:
summary(ephy, display = 10)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}

  \item{display}{
an integer for the number of rows of the posterior to display
}
}
\details{
Prints to console the number of posterior samples and the posterior 
distribution on the number of shifts, which is just the fraction of 
samples in the posterior having 0,1,2,...n shifts.
}
\value{
Returns (invisibly) a list with:
\item{posterior}{
	the number of posterior samples having 0,1,2,...,n shifts
}

}
\references{
\url{www.bamm-project.org}
}
\author{
Mike Grundler
}
 
\seealso{
}
\examples{
library(BAMMtools)
data(whales, events.whales)
ephy = getEventData(whales,events.whales,nsamples=100)
summary(ephy)
}

\keyword{ models }
