\name{plot.bammdata}
\alias{plot.bammdata}

\title{
Plot \code{BAMM}-estimated macroevolutionary rates on a phylogeny
}
\description{
\code{plot.bammdata} plots a phylogenetic tree from a \code{bammdata} object and
colors each branch by the estimated rate of speciation, extinction, or trait evolution. Rates are not assumed to be constant in time, and the function can plot continuously-varying rates along individual branches. 
}
\usage{
S3 method for class \code{bammdata}:		
plot(ephy, method = "phylogram", spex="s", 
    vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, multi = FALSE, 
    legend = FALSE, lwd = 1, cex = 1, pal = "set1", colorbreaks = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{method}{
a character string indicating the method for plotting the phylogenetic tree. 
\code{method = "phylogram"} (default) plots the phylogenetic tree using rectangular 
coordinates. \code{method = "polar"} plots the phylogenetic tree using polar coordinates.
}
  \item{spex}{
a character string indicating what type of macroevolutionary rates should be plotted. "s" (default) indicates
speciation rates, "e" indicates extinction rates, and any other character, e.g. "se", indicates net diversification
rates. Ignored if \code{ephy$type = "trait"}.  	
}
  \item{vtheta}{
a numeric indicating the angular separation (in degrees) of the first and last
terminal nodes. Ignored if \code{method = "phylogram"}.
}
  \item{rbf}{
a numeric indicating the length of the root branch as a fraction of total tree
height. Ignored if \code{method = "phylogram"}.
}
  \item{show}{
a logical indicating whether or not to plot the tree. Defaults to \code{TRUE}.
}
  \item{labels}{
a logical indicating whether or not to plot the tip labels. Defaults to \code{FALSE}.
}
  \item{multi}{
a logical indicating whether or not the device is configured for a multi-panel plot.
Defaults to \code{FALSE}.
}
  \item{legend}{
a logical indicating whether or not to plot a density function for interpreting the 
mapping of evolutionary rates to colors. Defaults to \code{FALSE}.
}
  \item{lwd}{
a numeric specifying the line width for branches.
}
  \item{cex}{
a numeric specifying the size of tip labels.
}
  \item{pal}{
a character string or vector of mode character that describes the color palette. See
Details for explanation of options.
}
  \item{colorbreaks}{
a numeric vector of percentiles delimiting the bins for mapping rates to colors. If \code{NULL} (default) bins
are calculated from the rates that are passed with the \code{bammdata} object.  	
}
  \item{\dots}{
further arguments passed to \code{par}.
}
}
\details{
To calculate rates, each branch of the phylogeny is discretized into a number of small segments, 
and the mean of the marginal posterior density of the rate of speciation/extinction or trait evolution is 
calculated for each such segment. Rates are mapped to colors such that cool colors represent 
slow rates and warm colors represent fast rates. When the tree is
plotted each of these small segments is plotted so that changes in rates through time and shifts in 
rates are visible as gradients of color.

Currently three options exist for specifying the color palette. First, any vector of three valid named colors
may be used. The first, second, and third colors are mapped to low, intermediate, and fast rates 
of evolution. This provides lots of flexibility for choosing a color scheme.
Second, the option \code{pal = "temperature"} uses the \code{rich.colors} function written by Arni Magnusson for the package gplots. Finally, if the package RColorBrewer is loaded, any of the named diverging color palettes documented 
under \code{brewer.pal} may also be used. Two default choices exist which make use of the first option. Specifying "set1"
uses a color ramp from green to red that we happen to like. This option is probably not suitable for people with certain forms of colorblindness. Specifying "set2" creates a color ramp from green to magenta that is hopefully better
for color deficient vision.

Specifying smaller values for \code{tau} will result in smoother-looking rate changes on the tree. Note that smaller values of \code{tau} require more computation.  

Internally \code{plot.bammdata} checks whether or not rates have been calculated by looking for a component named
"dtrates" in the \code{bammdata} object. If rates have not been calculated \code{plot.bammdata} calls \code{dtRates} with \code{tau = 0.01}. If the \code{colorbreaks} argument is \code{NULL} a map of rates to colors is also made by calling \code{assignColorBreaks} with \code{NCOLORS = 64}. A user supplied \code{colorbreaks} argument can be passed as well. This allows one to plot parts of a tree while preserving the map of rates to colors that was made using rates for the entire tree.

If \code{plot.bammdata} is called repeatedly with the same \code{bammdata} object computation can be reduced by first calling \code{dtRates} in the global environment.
}
\value{
Returns (invisibly) a list with three components. 
\item{coords}{
A matrix of plot coordinates. Rows correspond to branches. Columns 1-2 are starting (x,y) coordinates 
of each branch and columns 3-4 are ending (x,y) coordinates of each branch. If \code{method = "polar"} a 
fifth column gives the angle (in radians) of each branch.
}
\item{colorbreaks}{
A vector of percentiles used to group macroevolutionary rates into color bins.	
}
\item{colordens}{
A matrix of the kernel density estimates (column 2) of evolutionary rates (column 1) and the color (column 3) corresponding to each rate value.
}
}
\references{

}
\author{
Mike Grundler
}
\note{
If arguments to \code{par} are specified in \code{\dots} the function resets \code{par} variables before
returning unless \code{multi = TRUE}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dtRates}}, \code{\link{addBAMMshifts}}, \code{\link{assignColorBreaks}}, \code{\link{subtreeBAMM}}	
}
\examples{
	
library(BAMMtools)
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)

#the first call to plot.bammdata. No calculations or assignments of rates have been made
plot(ed, lwd = 3, spex = "s") #internally calls dtRates and assignColorBreaks

#now plot.bammdata no longer calls dtRates
ed <- dtRates(ed, tau = 0.01)
ed <- plot(ed, lwd = 3, spex = "s")

#you can plot subtrees while preserving the original rates to colors map by passing the colorbreaks object as an argument
sed <- subtreeBAMM(ed, node = 103)
plot(sed, lwd = 3, colorbreaks = ed$colorbreaks)
sed <- subtreeBAMM(ed, node = 140)
plot(sed, lwd = 3, colorbreaks = ed$colorbreaks)
#note how if we do not pass colorbreaks the map is no longer relative to the rest of the tree and the plot is quite
#distinct from the original
plot(sed, lwd = 3)

#if you want to change the value of tau and the rates to colors map for the entire tree
ed <- dtRates(ed, tau = 0.002)
#now you can re-plot the subtrees using this finer tau partition
sed <- subtreeBAMM(ed, node = 103)
sed <- dtRates(sed, 0.002)
plot(sed, lwd = 3, colorbreaks = ed$colorbreaks)
sed <- subtreeBAMM(ed, node = 140)
sed <- dtRates(sed, 0.002)
plot(sed, lwd = 3, colorbreaks = ed$colorbreaks)

#multi-panel plotting and adding shifts of specific posterior samples
par(mfrow=c(2,3),mar=c(5,1,1,1))
samples = sample(1:length(ed$eventData), 6)
for (i in 1:6) {
	ed <- dtRates(ed, 0.005, samples[i])
	plot(ed, multi=TRUE)
	addBAMMshifts(ed,method="phylogram",index=samples[i],multi=TRUE)	
}

#color options
ed <- dtRates(ed,0.01)
plot(ed, pal="set1")
plot(ed, pal="set2")
plot(ed, pal=c("midnightblue","lightskyblue","orangered"))
if (require(RColorBrewer)) {
	plot(ed,method="polar",pal="Spectral", lwd=3)
	plot(ed,method="polar",pal="RdYlBu", lwd=3)
}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
