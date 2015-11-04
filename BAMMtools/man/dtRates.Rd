\name{dtRates}
\alias{dtRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate macroevolutionary rate changes on a phylogeny from \code{BAMM} output
}
\description{
\code{dtRates} calculates the mean of the marginal posterior density 
of the rate of speciation/extinction or trait evolution for small segments
along each branch in a phylogeny.
}
\usage{
dtRates(ephy, tau, ism = NULL, tmat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{tau}{
a numeric that specifies the size (as a fraction of 
tree height) of the segments that each branch will be discretized into.
}
  \item{ism}{
an integer vector indexing which posterior samples to include
in the calculation.
}
  \item{tmat}{
  	a logical. If \code{TRUE} the matrix of branch segments is returned.
}
}
\details{
\code{dtRates} bins the phylogeny into windows of time and calculates average rates of speciation/extinction or phenotypic evolution along each segment of a branch within a window. The width of each window is determined by \code{tau}. \code{tau} is a fraction of the root to tip distance so a value of \code{tau = 0.01} bins the phylogeny into 100 time windows of equal width.
}

\value{
A \code{bammdata} object with a new component named "dtrates", which is a list with two or three components:	
\item{tau}{the parameter value of \code{tau} used in the calculation.}
\item{rates}{if \code{ephy$type = "trait"}: a numeric vector with the phenotypic rates of each segment on each branch. if \code{ephy$type = "diversification"}: a list with two components. The first component is a numeric vector of speciation rates. The second component is a numeric vector of extinction rates.}
\item{tmat}{a matrix of the starting and ending times of the segments on each branch. Only if \code{tmat = TRUE}.}
}
\references{
\url{bamm-project.org}
}
\author{
Mike Grundler
}
\note{
If there are zero length branches in the input tree \code{NA}s will result. 
}
 
\seealso{
\code{\link{plot.bammdata}}
}
\examples{
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)

#use all posterior samples
ed <- dtRates(ed, tau=0.01)

#use specified range of posterior samples
ed <- dtRates(ed, tau=0.01, ism=50:150)

}
 
\keyword{ graphics }
