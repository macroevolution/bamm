\name{assignColorBreaks}
\alias{assignColorBreaks}
\title{
	Map macroevolutionary rates to colors
}
\description{
	Maps macroevolutionary rates to a set of \code{NCOLORS}.
}
\usage{
	assignColorBreaks(NCOLORS, spex)
}
\arguments{
	\item{NCOLORS}{
		an integer number of colors to use for the mapping. Larger numbers do not 
		necessarily result in smoother looking color ramps. The default is 64 and 
		is probably sufficient for most purposes.
	}
	\item{spex}{
		a character string. "s" means that speciation rates are used to make the
		map, "e" means that extinction rates are used. Any other string, e.g. "se",
		means that diversification rates are used. Ignored for BAMM trait data.
	}
}
\details{
	Macroevolutionary rates are binned into \code{NCOLORS+1} percentiles and rates in each bin are
	mapped to a color determined by the \code{pal} argument in \code{plot.bammdata}.
}
\value{
	A "colorbreaks" object is assigned in the \code{.dtRatesEnv} environment 
}
\note{
	It does not make sense to call this function if \code{dtRates} has not been called and an error will
	result. The rates to colors map is relative to the set of posterior samples and nodes in \code{bammdata} object and the value of \code{tau} when \code{dtRates} was called.
}
\author{
	Mike Grundler
}
\seealso{
	\code{\link{plot.bammdata}}
}
\keyword{ graphics }