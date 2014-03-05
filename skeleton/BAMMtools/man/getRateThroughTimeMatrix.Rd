\name{getRateThroughTimeMatrix}
\alias{getRateThroughTimeMatrix}

\title{
Generate rate-through-time matrix from \code{bammdata} object
}
\description{
Computes a matrix of macroevolutionary rates at specified timepoints from a \code{bammdata} object. These rates can be used for plotting speciation rates (and other rates) through time. 
}
\usage{
getRateThroughTimeMatrix(ephy, start.time = NULL, 
	end.time = NULL, nslices = 100, node = NULL, nodetype = "include")
}

\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{start.time}{
the start time for the time sequence over which rates should be computed. If NULL, starts at time 0
}
  \item{end.time}{
the end time for the time sequence over which rates should be computed. If NULL, ends in the present
}
  \item{nslices}{
The number of time points at which to compute rate estimates (between \code{start.time} and \code{end.time})
}
  \item{node}{
Allows user to extract rate-through-time information for the subtree descended from a specific node. Alternatively, a specified subtree can be excluded from the rate matrix calculations.
}
  \item{nodetype}{
Two options: "include" and "exclude". If "include", computes rate matrix only for the descendants of subtree defined by \code{node}. If "exclude", computes rate matrix for all background lineages in tree after excluding the subtree defined by \code{node}. Ignored if \code{node = NULL}.
}

}
\details{
Computes evolutionary rates for each sample in the posterior included as part of the \code{bammdata} object. Rates are computed by draping an imaginary grid over the phylogeny, where the grid begins at \code{start.time} and ends at \code{end.time}, with \code{nslices} vertical lines through the phylogeny. The mean rate at each point in time (for a given sample from the posterior) is simply the mean rate at that time for all branches that are intersected by the grid (see the grid plot in the examples section).

This function is used by \link{plotRateThroughTime}, but the user can work directly with the \code{bamm-ratematrix} object for greater control in plotting rate-through-time trajectories for individual clades. See \code{examples} for an example of how this can be used to plot confidence intervals on a rate trajectory using shaded polygons.

The \code{node} options are particularly useful. If you have run \code{BAMM} on a large phylogeny, you can easily generate the rate-through-time data for a particular subtree by specifying the node number along with \code{nodetype = "include"}. Likewise, if you want to look at just the background rate - excluding some particular lineage - just specify \code{nodetype = "exclude"}. 

}
\value{
An object of class \code{bamm-ratematrix} with the following components:
\item{lambda}{ 
A \code{nsamples} x \code{nslices} matrix of speciation rates, where \code{nsamples} is the number of posterior samples in the \code{bammdata} object.
}
\item{mu}{ 
A \code{nsamples} x \code{nslices} matrix of extinction rates
}
\item{beta}{ 
A \code{nsamples} x \code{nslices} matrix of phenotypic rates (if applicable)
}
\item{times}{ 
A vector of timepoints where rates were computed
}
\item{times}{ 
A vector of timepoints where rates were computed (see example)
}
\item{type}{
Either "diversification" or "trait", depending on the input data. 
}

}
 
\author{
Dan Rabosky
}
 
\seealso{
\link[BAMMtools]{plotRateThroughTime}
}
\examples{
##NOT RUN##
# Plot a rate-through-time curve with 
# confidence intervals for the whale dataset:

# data(whales, events.whales)
# ed <- getEventData(whales, events.whales)

# rmat <- getRateThroughTimeMatrix(ed)

# plot.new()
# plot.window(xlim=c(0, 36), ylim=c(0, .7))

## Speciation quantiles: plot 90% CIs
# qq <- apply(rmat$lambda, 2, quantile, c(0.05, 0.5, 0.95))

# xv <- c(rmat$times, rev(rmat$times))
# yv <- c(qq[1,], rev(qq[3,]))

## Add the confidence polygon on rate distributions:
# polygon(xv, yv, col='gray80', border=FALSE)

## Add the median rate line:
# lines(rmat$times, qq[2,], lwd=3, col='red')

# Add axes
# axis(1, at=seq(-5, 35, by=5))
# axis(2, at=seq(-0.2, 1, by=0.2), las=1)

####### Now we will show the actual grid used for rate calculations:

# plot(whales, show.tip.label=FALSE)
# axisPhylo()

# mbt <- max(branching.times(whales))
## tvec <- mbt - rmat$times;
# tvec <- rmat$times;


# for (i in 1:length(tvec)){
# 	lines(x=c(tvec[i], tvec[i]), y=c(0, 90), lwd=0.7, col='gray70')
# }

## This shows the grid of time slices over the phylogeny
##END##
}
 
\keyword{ models}
