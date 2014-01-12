\name{getMeanTipRateStates}
\alias{getMeanTipRateStates}
\title{
Return mean tip rates
}
\description{
Return speciation, extinction, net diversification, or Brownian motion trait rates for all tips in the phylogeny.
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
boolean: should the returned vector be named with taxon names, defaults to \code{TRUE}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
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
library(BAMMtools)
data(cetaceans, events.cetaceans)
ephy <- getEventData(cetaceans, events.cetaceans, burnin=0.25)

#return a vector of net diversification rates at the tips of the phylogeny.
getMeanTipRateStates(ephy, ratetype='netdiv')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
