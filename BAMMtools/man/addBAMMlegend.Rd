\name{addBAMMlegend}
\alias{addBAMMlegend}
\title{
Add a color legend to a phylo-rate plot
}
\description{
Add a legend to a phylorate plot, with greater manual control.
}
\usage{
addBAMMlegend(x, corners = c(0,1,0,10), side = 3, nTicks = 2, 
	direction = "auto", axisOffset = 0.4, ...)
}
\arguments{
	\item{x}{
		a plot.bammdata object.
	}
	\item{corners}{
		coordinates for the corners of the bar legend \code{c(xmin, xmax, ymin, ymax)}.	
	}
	\item{side}{
		side for tick marks, see \code{\link{axis}} documentation.
	}
	\item{nTicks}{
		number of tick marks, besides of min and max.
	}
	\item{direction}{
		direction of color ramp. If \code{auto}, then direction is automatically inferred,
		otherwise can be specified as \code{horizontal} or \code{vertical}.
	}
	\item{axisOffset}{
		distance from color bar for labels.
	}
	\item{...}{
		additional parameters to be passed to \code{\link{axis}}.
	}
}

\value{
	Legend is plotted. Nothing is returned.
}
\author{
	Pascal Title
}
\seealso{
	Requires an object created with \code{\link{plot.bammdata}}. 
}
\examples{
data(whales, events.whales)
ephy <- getEventData(whales, events.whales, burnin=0.25)

x <- plot.bammdata(ephy, lwd=2)
axis(1); axis(2)
addBAMMlegend(x, corners = c(-1, -0.5, 40, 80), nTicks = 4, side = 4, las = 1)
addBAMMlegend(x, corners = c(5, 20, 60, 61), nTicks = 4, side = 3, las = 1, cex.axis = 0.7)

#addBAMMlegend also automatically detects the use of color.interval
data(primates, events.primates)
ephy <- getEventData(primates, events.primates, burnin=0.25, type = 'trait')

x <- plot.bammdata(ephy, breaksmethod='linear', color.interval=c(NA, 0.12), lwd=2)
addBAMMlegend(x, corners = c(0, 30, 200, 205), nTicks = 1, side = 3, axisOffset = 1, las = 1)
}
 

