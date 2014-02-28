\name{plot.credibleshiftset}
\alias{plot.credibleshiftset}

\title{
	Plot core rate shift configurations on a phylogeny
}
\usage{
S3 method for class \code{credibleshiftset}
plot(cset, plotmax=9, method="phylogram", pal="RdYlBu", spex="s", add.freq.text=TRUE,
use.plot.bammdata=TRUE, border=TRUE, legend=FALSE, send2pdf=FALSE, ...)
}
\arguments{
\item{cset}{
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
	
}
\value{
	
}
\references{
\link{bamm-project.org}	
}
\author{
Mike Grundler
}
\seealso{
\code{\link{credibleShiftSet}}, \code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammdata}}
}
\examples{
data(events.whales)
data(whales)
ed <- getEventData(whales, events.whales)
cset <- credibleShiftSet(ed)
plot(cset)
}
