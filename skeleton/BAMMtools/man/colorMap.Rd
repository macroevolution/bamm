\name{colorMap}
\alias{colorMap}
\title{
Map rates to colors
}
\description{
Maps rates of speciation or trait evolution to color space.
}
\usage{
colorMap(x, pal, NCOLORS)
}
\arguments{
	\item{x}{
a numeric vector of rates of speciation or trait evolution.		
}
	\item{pal}{
a character string or vector of mode character specifying the palette. See 
documentation of \code{plot.bammdata} for further details.		
}
	\item{NCOLORS}{
an integer number of colors used to create the color gradient.		
}
}
\details{
Rates of speciation or trait evolution are mapped to colors by first creating
a color gradient, \code{colpalette}, of \code{NCOLORS} colors. The vector \code{x} 
of rates is then binned into a vector, \code{breaks}, of \code{NCOLORS+1} percentiles. 
Values of \code{x} in the interval [\code{breaks[i-1]},\code{breaks[i]}) are assigned 
to color \code{colpalette[i-1]}.
}
\value{
A vector of \code{length(x)} colors.
}
\references{
	
}
\author{
	Dan Rabosky, Mike Grundler
}
\note{
	
}
\seealso{
\code{\link{plot.bammdata}}, \code{\link{dtRates}}.
}
\examples{
	
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }