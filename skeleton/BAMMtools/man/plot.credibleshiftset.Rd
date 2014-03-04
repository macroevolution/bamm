\name{plot.credibleshiftset}
\alias{plot.credibleshiftset}

\title{
	Plot credible set of rate shift configurations from BAMM analysis
}
\description{
	Plots the credible set of rate shift configurations from a BAMM analysis on a phylogeny
}
\usage{
\method{plot}{credibleshiftset}(x, plotmax=9, method="phylogram", pal="RdYlBu", spex="s", add.freq.text=TRUE,
use.plot.bammdata=TRUE, border=TRUE, legend=FALSE, send2pdf=FALSE, ...)
}
\arguments{
\item{x}{
an object of class \code{credibleshiftset}.	
}
\item{plotmax}{
an integer number of plots to display. 
}
\item{method}{
a coordinate method to use for plotting. Options are "phylogram" or "polar".
}
\item{pal}{
a color palette to use with \code{plot.bammdata}.
}
\item{spex}{
a character string indicating what type of macroevolutionary rates should be plotted. "s" (default) indicates speciation rates, "e" indicates extinction rates, and any other character, e.g. "se", indicates net diversification rates. Ignored if ephy$type = "trait".
}
\item{add.freq.text}{
a logical indicating whether to add the posterior frequency of each shift configuration to the plotting region.	
}
\item{use.plot.bammdata}{
a logical indicating whether to use \code{plot.bammdata} (\code{TRUE}) or \code{plot.phylo} (\code{FALSE}).	
}
\item{border}{
a logical indicating whether to frame the plotting region.	
}
\item{legend}{
a logical indicating whether to plot a legend.
}
\item{send2pdf}{
a logical indicating whether to print the figure to a PDF file.
}
\item{...}{
further arguments to pass to \code{plot.bammdata}.
}
}
\details{
This produces phylorate plots for the \code{plotmax} most-probable shift configurations sampled with \code{BAMM}. Shift configurations are plotted in a single graphics window. The posterior probability (frequency) of each rate shift configuration in the posterior is shown (omitted with argument \code{add.freq.text = FALSE}).

Points are added to the branches subtending the nodes of each rate configuration. The size of the point is proportional to the marginal probability that a shift occurs on a specific branch. If the instantaneous rate at a shift's origin represents an initial increase above the ancestral instantaneous rate the point is red. If the instantaneous rate at a shift's origin represents an initial decrease below the ancestral instantaneous rate the point is blue.  	
}

\references{
\url{bamm-project.org}	
}
\author{
Mike Grundler
}
\seealso{
\code{\link{credibleShiftSet}}, \code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammdata}}, \code{\link{plot.bammshifts}}
}
\examples{
data(events.whales)
data(whales)
ed <- getEventData(whales, events.whales)

# Get branch prior on number of shifts:
data(prior.whales)
bp <- getBranchShiftPriors(whales, prior.whales)
cset <- credibleShiftSet(ed, threshold = bp)
plot(cset)
}
