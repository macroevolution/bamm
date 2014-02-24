\name{credibleShiftSet}
\alias{credibleShiftSet}
\title{
Credible set of macroevolutionary rate shift configurations from BAMM results
}
\description{
Computes the XX\% credible set of macroevolutionary rate shift configurations from a \code{bammdata} object. These results can be analyzed further and/or plotted.

}
\usage{
credibleShiftSet(ephy, set.limit = 0.95, threshold = 0.01)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{set.limit}{
the desired limit to the credible set. A value of 0.95 will return the 95\% credible set of shift configurations.
}
  \item{threshold}{
Marginal probability threshold for including nodes in the enumeration of distinct shift configurations. Recommended value = 0.01. See details.
}
}
\details{
Computes the XX\% credible set of diversification shift configurations sampled using BAMM. This is analogous to a credible set of phylogenetic tree topologies from a Bayesian phylogenetic analysis. 

To understand how this calculation is performed, one must first distinguish between "core" and "transient" rate shifts. Because the BAMM algorithm is always proposing (and deleting) shifts during simulation of the posterior, it is usually the case that - for any sample - there is one or more "transient" shifts that has been proposed and accepted but which is not really improving the ability to model the data. Hence, such transient shifts will have extremely low marginal probabilities overall. Core shifts are those that do contribute to the ability to model the data and which are favored in the posterior as they improve the posterior probability of the data. 

If we compute distinct shift configurations with every sampled shift (including those shifts with very low marginal probabilities), the number of distinct shift configurations will be overwhelmingly high. However, most of these configurations are functionally identical because they contain the same core shifts and differ only in transient shifts that have very low marginal probabilities. 

By specifying a \code{threshold} argument, the user can choose to ignore shifts with low probabilities during enumeration of distinct shift configurations. We suggest using this option. Rate shifts with marginal probabilities less than \code{threshold} are simply ignored. 
 
}
\value{
A class \code{credibleshiftset} object with the following components:
  \item{shiftnodes}{
A list of the "core" rate shifts (marginal probability > threshold) that occurred in each distinct shift configuration in the credible set. The i'th vector from this list gives the core shift nodes for the i'th shift configuration. They are sorted by frequency, so \code{x\$shiftnodes[[1]]} gives the shift nodes that occurred together in the shift configuration with the highest posterior probability.
}
  \item{indices}{
A list of vectors containing the indices of samples in the \code{bammdata} object that are assigned to a given shift configuration. All are sorted by frequency.
}
  \item{frequency}{
A vector of frequencies of shift configurations, including those that account for \code{set.limit} (typically, 0.95 or 0.99) of the probability of the data.
}
  \item{cumulative}{
Like \code{frequency}, but contains the cumulative frequencies
}
  \item{number.distinct}{
Number of distinct shift configurations in the credible set
}
  \item{set.limit}{
which credible set is this (0.9, 0.95, etc)?
}
}
\references{
\code{bamm-project.org}
}
\author{
Dan Rabosky
}

\seealso{
\code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammshifts}}, \code{\link{summary.credibleshiftset}}	
} 
 
 
\examples{
data(events.whales)
data(whales)
ed <- getEventData(whales, events.whales)
cset <- credibleShiftSet(ed)
summary(cset) 

}
 
\keyword{ models }
