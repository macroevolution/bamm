\name{distinctShiftConfigurations}
\alias{distinctShiftConfigurations}

\title{
Identify distinct rate shift configurations
}
\description{
Identify topologically distinct rate shift configurations that were sampled with BAMM, and assign each sample in the posterior to one of the distinct shift configurations.	
}
\usage{
distinctShiftConfigurations(ephy, threshhold)
}
\arguments{
	\item{ephy}{
		an object of class \code{bammdata}.
	}
	\item{threshold}{
		Either an object of class \code{branchprior} or a numeric indicating the minimum frequency at which a shift has to be sampled to be included in the calculation.
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
The marginal probability threshold for nodes (ignores all nodes less than this value). Alternatively, an object of class \code{branchprior} containing this information.	
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
	\code{\link{plot.bammshifts}}, \code{\link{credibleShiftSet}}
}
\examples{
library(BAMMtools)
data(whales, events.whales, prior.whales)

# Get the prior distribution on the number of shifts per branch:
bp <- getBranchShiftPriors(whales, prior.whales)


ed <- getEventData(whales, events.whales, burnin=0.25)

sc <- distinctShiftConfigurations(ed, threshold=bp)

plot.bammshifts(sc, ed, rank=1)

}



