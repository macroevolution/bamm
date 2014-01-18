\name{maximumShiftCredibility}
\alias{maximumShiftCredibility}
\title{
Compute maximum shift credibility configuration
}
\description{
Finds the rate shift configuration with the highest marginal probability. It is one estimate of the "best" rate shift configuration, considering only those shift configurations that were actually sampled using \code{BAMM}'s reversible jump MCMC simulator. This is analogous to the "maximum clade credibility tree" from a Bayesian phylogenetic analysis. 

}
\usage{
maximumShiftCredibility(ephy, maximize = "product")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
an object of class \code{BAMMdata}
}
  \item{maximize}{
Maximize the marginal probability of the product or sum of branch-specific shifts (see details)
}
}
\details{
This is one point estimate of the overall "best" rate shift configuration. Following an MCMC simulation, the marginal shift probabilities on each individual branch are computed using \link[BAMMtools]{marginalShiftProbsTree}. The shift configuration that maximizes the product (or sum, if specified) of these marginal branch-specific shift probabilities is the \emph{maximum shift credibility configuration}.  
}
\value{
A list with the following components:
\item{bestconfigs}{ 
A vector of the index values of MCMC samples with shift configurations equal to the maximum. Usually, more than one state sampled during the MCMC simulation will have an identical (maximized) marginal probability. All samples given in this vector will have an identical shift configuration.
}

\item{scores}{The optimality score (product or sum of marginal shift probabilities) for all sampled shift configurations in the \code{BAMMdata} object
}

\item{optimalityType}{Whether the product or sum of marginal shift probabilities was used to compute the maximum shift credibility configuration
} 

\item{sampleindex}{A representative sample that is equal to the maximum shift credibility configuration (e.g., this can be plotted with \link[BAMMtools]{addBAMMshifts})
}
}

\references{

}

\author{
Dan Rabosky
}

\seealso{
\link[BAMMtools]{marginalShiftProbsTree}, \link[BAMMtools]{addBAMMshifts}, \link{cumulativeShiftProbsTree}
}
\examples{
library(BAMMtools)
data(whales)
data(events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)
best_config <- maximumShiftCredibility(ed)
plot.bammdata(ed)
addBAMMshifts(ed, method='phylogram', index=best_config$sampleindex)
}

\keyword{ manip }
\keyword{graphics}
