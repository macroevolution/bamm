\name{distinctShiftConfigurations}
\alias{distinctShiftConfigurations}

\title{
	Sample distinct rate shift configurations from the posterior
}
\description{
	
}
\usage{
distinctRateShiftConfigurations(ephy, threshhold=0.05)
}
\arguments{
	\item{ephy}{
		an object of class \code{bammdata}.
	}
	\item{threshold}{
		a numeric indicating the minimum frequency at which a shift has to be sampled
		to be included in the calculation.
	}
}
\details{
	
}
\value{
an object of class \code{bammshifts}. This is a list with the following components
\item{marg.probs}{
A list of the marginal probability of a shift occurring at each node of the phylogeny for each
distinct rate configuration.	
}
\item{shifts}{
A list of the set of shift nodes for each distinct rate configuration.	
}
\item{samplesets}{
A list of sample indices that reduce to each of the unique shift sets.	
}
\item{frequency}{
A vector of frequencies of each distinct shift configuration	.
}
\item{threshold}{
The marginal probability threshold for nodes (ignores all nodes less than this value).	
}	
Results are sorted by frequency:
$frequency[1] gives the most common shift configuration sampled.
$shifts[[1]] gives the corresponding node indices for that configuration.
$samplesets[[1]] gives the indices of samples with this configuration.
}
\references{
	
}
\author{
	Dan Rabosky
}
\seealso{
	\code{\link{plot.bammshifts}}
}
\examples{
library(BAMMtools)
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)

sc <- distinctShiftConfigurations(ed, threshold=0.05)
plot(sc, ed)
}