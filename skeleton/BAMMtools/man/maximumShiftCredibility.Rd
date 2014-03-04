\name{maximumShiftCredibility}
\alias{maximumShiftCredibility}
\title{
Estimate maximum shift credibility configuration
}
\description{
This is one estimate of the "best" rate shift configuration, considering only those shift configurations that were actually sampled using \code{BAMM}'s reversible jump MCMC simulator. This is analogous to the "maximum clade credibility tree" from a Bayesian phylogenetic analysis. It is not necessarily the same as the shift configuration with the maximum a posteriori probability.

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
Maximize the marginal probability of the product or sum of branch-specific shifts
}
}
\details{
This is one point estimate of the overall "best" rate shift configuration. Following an MCMC simulation, the marginal shift probabilities on each individual branch are computed using \link[BAMMtools]{marginalShiftProbsTree}. The shift configuration that maximizes the product (or sum, if specified) of these marginal branch-specific shift probabilities is the \emph{maximum shift credibility configuration}.  

This option is only recommended if you have no clear "winner" in your credible set of shift configurations (see \link[BAMMtools]{credibleShiftSet}). If you have a number of largely-equiprobable shift configurations in your 95\% credible set, you may wish to try this function as an alternative for identifying a single best shift configuration. Otherwise, it is recommended that you present the shift configuration with the maximum a posteriori probability (see \link[BAMMtools]{getBestShiftConfiguration}).

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

\author{
Dan Rabosky
}

\seealso{
\link[BAMMtools]{marginalShiftProbsTree}, \link[BAMMtools]{addBAMMshifts}, \link{cumulativeShiftProbsTree}, \link[BAMMtools]{credibleShiftSet}, \link[BAMMtools]{getBestShiftConfiguration}
}
\examples{
data(whales)
data(events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)
best_config <- maximumShiftCredibility(ed)
plot.bammdata(ed)
addBAMMshifts(ed, method='phylogram', index=best_config$sampleindex)
}

\keyword{ manip }
\keyword{graphics}
