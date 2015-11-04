\name{getTipRates}
\alias{getTipRates}
\title{
Compute tip-specific macroevolutionary rates from \code{bammdata} object
}
\description{
Return speciation, extinction, net diversification, or Brownian motion trait rates for all species in the phylogeny from \code{BAMM} output.
}
\usage{
getTipRates(ephy, returnNetDiv = FALSE, statistic = 'mean')
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{returnNetDiv}{
logical. If \code{TRUE}, then net diversification rates are returned, if \code{FALSE}, then both speciation and extinction rates are returned. If \code{ephy} is of type \code{trait}, then this is ignored.
}
  \item{statistic}{
Determines how the average tip rates should be calculated. Can be either \code{mean} or \code{median}.
}
}

\value{
Returns a list with the following elements:
If \code{ephy} type is 'diversification':
\item{lambda}{
	A matrix of tip speciation rates with species as rows, and posterior samples as columns.
}
\item{mu}{
	A matrix of tip extinction rates with species as rows, and posterior samples as columns.
}
\item{lambda.avg}{
	A vector of average tip speciation rates, averaged with mean or median, depending on selected option for \code{statistic}. The vector is named with species names. 
}
\item{mu.avg}{
	A vector of average tip extinction rates, averaged with mean or median, depending on selected option for \code{statistic}. The vector is named with species names. 
}

If \code{ephy} type is 'diversification' and \code{returnNetDiv = TRUE}:
\item{netdiv}{
	A matrix of tip net diversification rates with species as rows, and posterior samples as columns.
}
\item{netdiv.avg}{
	A vector of average tip net diversification rates, averaged with mean or median, depending on selected option for \code{statistic}. The vector is named with species names. 
}

If \code{ephy} type is 'trait':
\item{beta}{
	A matrix of tip phenotypic rates with species as rows, and posterior samples as columns.
}
\item{beta.avg}{
	A vector of average tip phenotypic rates, averaged with mean or median, depending on selected option for \code{statistic}. The vector is named with species names. 
}


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

#return a vector of average species-specific speciation rates.
meanlam <- getTipRates(ephy, returnNetDiv = FALSE, statistic = 'mean')$lambda.avg
meanlam

#return a vector of median species-specific net diversification rates.
ndr <- getTipRates(ephy, returnNetDiv = TRUE, statistic = 'median')$netdiv.avg


#Return mean species-specific speciation rates from all posterior 
#	samples in the \code{bamm-data} object.
lam <- getTipRates(ephy, returnNetDiv = FALSE, statistic = 'mean')$lambda
rowMeans(lam)


}
 
\keyword{ models }
