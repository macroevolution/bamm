\name{summary.bammdata}
\alias{summary.bammdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary
%%  ~~function to do ... ~~
}
\description{
Summarizes the number of shifts and their locations in the 
maximum shift credibility (MSC) tree and the posterior 
distribution on the number of shifts.
}
\usage{
S3 method for class \code{bammdata}:
summary(ephy, opt = "product", display = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{opt}{
a character string indicating the optimality criterion. May
be "product" or "sum"
}
  \item{display}{
an integer for the number of rows of the posterior to display
}
}
\details{
Prints to console the number of posterior samples, the number of 
shifts in the MSC tree, the locations of those shifts, the optimality
criterion used for calculating the MSC tree, and the posterior 
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
	tree
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mike Grundler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{maximumShiftCredibilityTree}}
}
\examples{
library(BAMMtools)
data(cetaceans, events.cetaceans)
ephy = getEventData(cetaceans,events.cetaceans,nsamples=100)
summary(ephy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
