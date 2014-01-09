\name{plot.bammdata}
\alias{plot.bammdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting dynamic rates on a phylogeny
}
\description{
\code{plot.bammdata} plots the phylogenetic tree in the \code{bammdata} object.
In doing so branches are colored according to estimated dynamic rates of speciation 
or trait evolution.	 
}
\usage{
S3 method for class \code{bammdata}:		
plot(ephy, method = "phylogram", tau = 0.01, index = NULL, 
    vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, multi = FALSE, 
    hrates = FALSE, lwd = 1, cex = 1, ncolors = 64, pal = "Spectral", ...)
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
  \item{tau}{
a numeric indicating the step size used in discretizing branches into segments. 
These segments are used in calculating rate changes along branches. Should be
given as a fraction of the total tree height, e.g. 0.01 (default). 
}
  \item{index}{
an integer vector indexing which posterior samples to include in the calculation
of rates. If \code{NULL} (default) all posterior samples in \code{ephy} are used.
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
  \item{hrates}{
a logical indicating whether or not to plot a density function for interpreting the 
mapping of evolutionary rates to colors. Defaults to \code{FALSE}.
}
  \item{lwd}{
a numeric specifying the line width for branches.
}
  \item{cex}{
a numeric specifying the character expansion for tip labels.
}
  \item{ncolors}{
an integer specifying the number of colors to be used for the mapping of evolutionary
rates to colors.
}
  \item{pal}{
a character string or vector of mode character that describes the color palette.
}
  \item{\dots}{
further arguments passed to \code{par}.
}
}
\details{
Internally \code{plot.bammdata} looks for a list element in \code{ephy} called "dtrates". 
This is simply a vector holding the rates of each segment on each branch. If "dtrates" is not 
found \code{plot.bammdata} calls \code{dtRates}, which performs the calculations of rates.
If you plan to repeatedly call \code{plot.bammdata} on a large tree you can reduce computation
load by calling \code{dtRates} with \code{ephy} before passing \code{ephy} to \code{plot.bammdata}.

Currently three options exist for specifying the color palette. In each case cool colors represent 
slow rates of speciation or trait evolution and warm colors represent fast rates. The default setting 
\code{pal = "temperature"} uses the \code{rich.colors} function from the package ggplots. If the package 
RColorBrewer is available then any of the named diverging color palettes documented under \code{brewer.pal} 
may also be used. Finally a vector of three valid named colors (e.g. c("blue","green","red")) may also be 
used. The first, second, and third colors are mapped to low, intermediate, and fast rates of evolution. 
}
\value{
Returns (invisibly) a matrix. Rows correspond to branches. Columns 1-2 are starting (x,y) coordinates 
of each branch and columns 3-4 are ending (x,y) coordinates of each branch. If \code{method = "polar"} a 
fifth column gives the angle (in radians) of each branch.
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
\code{\link{dtRates}}.	
}
\examples{
	
library(BAMMtools)
data(cetaceans, events.cetaceans)
ephy = getEventData(cetaceans, events.cetaceans, burnin=0.25)

#adding shifts to tree for specific posterior samples
plot(ephy, method="polar", index=5)
addBAMMshifts(ephy, "polar", index=5)

#multi-panel plotting and adding shifts
par(mfrow=c(2,3),mar=c(5,1,1,1))
samples = sample(1:length(ephy$eventData), 6)
for (i in 1:6) {
	plot(ephy, tau=0.001, index=samples[i], multi=TRUE)
	addBAMMshifts(ephy,method="phylogram",index=samples[i],multi=TRUE)	
}

#default behavior will call dtRates on each invocation
plot(ephy)
plot(ephy, method="polar")

#now plot.bammdata does not call dtRates
ephy = dtRates(ephy)
plot(ephy)
plot(ephy,method="polar")

#color and line width options
if (require(RColorBrewer)) {
	plot(ephy,method="polar",pal="Spectral", lwd=3)
	plot(ephy,method="polar",pal="RdYlBu", lwd=3)
}	
plot(ephy,method="polar",pal=c("blue","green","red"), lwd=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
