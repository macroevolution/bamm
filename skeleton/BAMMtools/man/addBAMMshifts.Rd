\name{addBAMMshifts}
\alias{addBAMMshifts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add BAMM shifts to a plotted tree
}
\description{
Adds symbols to a plotted tree to mark the location(s) where there is a shift in the
macroevolutionary dynamics of diversification or trait evolution.
}
\usage{
addBAMMshifts(ephy, method, index, cex = 1, pch = 21, col = 1, bg = 2, multi = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{method}{
a character string indicating the method used in plotting. May be "polar" or "phylogram".
}
  \item{index}{
an integer indicating which posterior sample to use for adding shifts to the plotted tree.
}
  \item{cex}{
a numeric indicating the character expansion of the plotted points.
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
  \item{multi}{
a logical indicating whether the device is configured for multi-panel plotting.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mike Grundler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getShiftNodesFromIndex}}, \code{\link{plot.bammdata}}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
