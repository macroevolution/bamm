\name{dtRates}
\alias{dtRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate macroevolutionary rate changes on a phylogeny from \code{BAMM} output
}
\description{
\code{dtRates} calculates the mean of the marginal posterior density 
of the rate of speciation or trait evolution for small segments
along each branch in a phylogeny.
}
\usage{
dtRates(ephy, tau, ism = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{tau}{
a numeric that specifies the approximate size (as a fraction of 
tree height) of the segments that each branch will be discretized into.  
}
  \item{ism}{
an integer vector indexing which posterior samples to include
in the calculation.
}
}
\details{
Rates of speciation and trait evolution are modeled by an exponential change 
function \eqn{R(t)=R(0)*exp(b*t)}. \eqn{R(0)} is the initial rate and \eqn{b}
determines the rate of increase or decrease of the rate \eqn{R(t)} through time.  
The average rate of speciation or trait evolution along a branch or segment of a 
branch with starting and ending times \eqn{t1} and \eqn{t2}, respectively, is: 

1/\eqn{(t2-t1)} * integral\eqn{(R(t), t1, t2) =} 1/\eqn{(t2-t1) * R(0)/b * (exp(b*t2) - exp(b*t1))}

\code{dtRates} will discretize each branch into a number of segments determined by
the parameter \code{tau}. Specifically, if \eqn{x} equals the branch length (as a 
fraction of tree height) and \eqn{n} equals the number of segments on that branch, 
\eqn{n} is determined by:

if (\eqn{x/}\code{tau} mod 1 == 0) \eqn{n = x/}\code{tau} + 1 else \eqn{n = ceiling(x/}\code{tau}\eqn{)}

This is why \code{tau} specifies the approximate size of the segments, since they will
differ from \code{tau} by a small amount across branches.
}

\value{
A \code{bammdata} object with a new component named "dtrates", which is a list with two components: 
\item{tau}{the parameter value of \code{tau} used in the calculation.}
\item{rates}{a numeric vector with the rates of each segment on each branch.}
}
\references{

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
library(BAMMtools)
data(whales, events.whales)
ephy = getEventData(whales, events.whales, burnin=0.25)

#use all posterior samples
ephy = dtRates(ephy, tau=0.01)

#use specified range of posterior samples
ephy = dtRates(ephy, tau=0.01, ism=50:150)
}
 
\keyword{ graphics }
