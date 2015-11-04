\name{getBestShiftConfiguration}
\alias{getBestShiftConfiguration}

\title{
Get the best (sampled) rate shift configuration from a BAMM analysis
}

\description{
Get the rate shift configuration with the maximum a posteriori probability, e.g., the shift configuration that was sampled most frequently with \code{BAMM}.
}
\usage{
getBestShiftConfiguration(x, expectedNumberOfShifts, threshold  = 5)
}

\arguments{
  \item{x}{
either a \code{bammdata} object or a \code{credibleshiftset} object
}
\item{expectedNumberOfShifts}{
The expected number of shifts under the prior
}
\item{threshold}{
The marginal posterior-to-prior odds ratio used as a cutoff for distinguishing between "core" and "non-core" rate shifts.
}
}
\details{
This function estimates the rate shift configuration with the highest maximum a posteriori (MAP) probability. It returns a \code{bammdata} object with a single sample. This can be plotted with \link[BAMMtools]{plot.bammdata}, and individual rate shifts can then be added with \link[BAMMtools]{addBAMMshifts}. 

The parameters of this object are averaged over all samples in the posterior that were assignable to the MAP shift configuration. All non-core shifts have been excluded, such that the only shift information contained in the object is from the "significant" rate shifts, as determined by the relevant marginal posterior-to-prior odds ratio \code{threshold}.

You can extract the same information from the credible set of shift configurations. See \link[BAMMtools]{credibleShiftSet} for more information.

}
\value{
A class \code{bammdata} object with a single sample, corresponding to the diversification rate shift configuration with the maximum a posteriori probability. See \link[BAMMtools]{getEventData} for details.
}
\author{
Dan Rabosky
}

\seealso{
\link{getEventData},  \link{credibleShiftSet}, \link{plot.credibleshiftset}, \link{plot.bammdata}
}
\examples{
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.1)

# Get prior distribution on shifts-per-branch:
bp <- getBranchShiftPriors(whales, expectedNumberOfShifts = 1)

# Pass the event data object in to the function:
best <- getBestShiftConfiguration(ed, expectedNumberOfShifts = 1, threshold = 5)
plot.bammdata(best, lwd=2)
addBAMMshifts(best, cex=2)

# Now we can also work with the credible shift set:
css <- credibleShiftSet(ed, expectedNumberOfShifts = 1, threshold = 5) 

summary(css)

# examine model-averaged shifts from MAP configuration-
# This gives us parameters, times, and associated nodes
#	of each evolutionary rate regime 
#   (note that one of them corresponds to the root)
css$eventData[[1]];

# Get bammdata representation of MAP configuration:
best <- getBestShiftConfiguration(css, expectedNumberOfShifts = 1, threshold = 5);

plot.bammdata(best)
addBAMMshifts(best)

 
 
}
\keyword{ models }
