\name{addBAMMshifts}
\alias{addBAMMshifts}

\title{
Add BAMM-inferred rate shifts to a phylogeny plot
}
\description{
Adds symbols to a plotted tree to mark the location(s) where there is a shift in the
macroevolutionary dynamics of diversification or trait evolution.
}
\usage{
addBAMMshifts(ephy, index = 1, method = 'phylogram', 
	cex = 1, pch = 21, col = 1, bg = 2, msp = NULL, 
		shiftnodes = NULL, par.reset = TRUE)
}

\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{method}{
a character string indicating the method used in plotting. Must be "polar" or "phylogram".
}
  \item{index}{
an integer indicating which posterior sample to use for adding shifts to the plotted tree.
}
  \item{cex}{
a numeric indicating the character expansion ("size") of the plotted points.
}
  \item{pch}{
an integer indicating the choice of plotting symbol.
}
  \item{col}{
an integer or character string indicating the border color of the plotting symbol.
}
  \item{bg}{
an integer or character string indicating the background color of the plotting symbol.
}
  \item{msp}{
if not \code{NULL}, an object of class \code{phylo} where each branch length is equal to 
the marginal probability of a shift occurring on that branch. Plotted points corresponding
to shifts will be sized by these probabilities.  	
}
  \item{shiftnodes}{
an optional vector of node numbers indicating the locations of shifts to plot. 	
}
  \item{par.reset}{
a logical indicating whether to reset the graphical parameters before exiting.
}
}
\details{

Any given sample from the posterior distribution sampled using \code{BAMM} contains a potentially unique configuration of rate shifts and associated parameters. There is no single "best" rate shift, but rather a set of shift configurations (and associated parameters) - along with their relative probabilities - sampled with MCMC. This function enables the user to plot the locations of shifts sampled with \code{BAMM} for a given sample from the posterior. 

If the \code{bammdata} object contains just a single sample, these shifts will be plotted regardless of the value of \code{index}.


}

\note{
If a \code{shiftnodes} argument is passed care should be taken to ensure that the nodes are in the same order
as in the event data for the sample index. 	
}

\author{
Mike Grundler
}

\seealso{
\code{\link{getShiftNodesFromIndex}}, \code{\link{plot.bammdata}}
}
\examples{

data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)

#adding shifts to tree for specific posterior samples
plot(ed, method="polar")
addBAMMshifts(ed, index=5, "polar")

#multi-panel plotting and adding shifts
par(mfrow=c(2,3),mar=c(5,1,1,1))
samples = sample(1:length(ed$eventData), 6)
for (i in 1:6) {
  sed <- subsetEventData(ed, samples[i])
  plot(sed, par.reset=FALSE)
  addBAMMshifts(sed,index=1,method="phylogram",par.reset=FALSE)	
}
}

\keyword{ graphics }
