\name{ratesHistogram}
\alias{ratesHistogram}

\title{
	Histogram of BAMM rate frequencies
}
\description{
	Plots a histogram of the frequency of rate values across the phylogeny.
}
\usage{
ratesHistogram(phylorates, plotBrks = TRUE, xlab = 'speciation rate', ylab = 'density', lwd = 0.2, lty = 1, brksCol = 'black', ...) 
}
\arguments{
	\item{phylorates}{
		a saved \code{\link{plot.bammdata}} object
	}
	\item{plotBrks}{
		boolean, should breaks be plotted over the histogram
	}
	\item{xlab}{
		x-axis label
	}
	\item{ylab}{
		y-axis label
	}
	\item{lwd}{
		line width for breaks
	}
	\item{lty}{
		line style for breaks
	}
	\item{brksCol}{
		color of breaks lines
	}
	\item{...}{
		additional arguments passed on to \code{mtext} for axis labels
	}
}
\details{
	With this function, a histogram is plotted that shows the frequency of rates present in the dataset. The color scheme plotted is taken from the saved \code{plot.bammdata} object that is the main input. Therefore, the mapping of colors to rates in the histogram corresponds exactly to what is plotted in the phylorate plot. 
	If \code{plotBrks = TRUE}, then the color breaks used for the phylorates plot are shown. 

	This function can be a useful tool for exploring different \code{plot.bammdata} options. Please see \url{http://bamm-project.org/colorbreaks.html} on the bamm-project website for more information on the utility of this function. 
}

\references{
	\url{bamm-project.org}
}
\author{
	Pascal Title
}

\seealso{
	\code{\link{plot.bammdata}}
}

\examples{

	library(BAMMtools)

	#load data
	data(primates, events.primates)
	ed <- getEventData(primates, events.primates, burnin=0.25, type = 'trait')

	#create phylorate plot with the jenks breaks method to generate output
	phylorates <- plot.bammdata(ed, breaksmethod='jenks', show = FALSE)

	ratesHistogram(phylorates, plotBrks = TRUE, xlab = 'trait rates')

}




