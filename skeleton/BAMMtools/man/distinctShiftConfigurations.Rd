\name{distinctShiftConfigurations}
\alias{distinctShiftConfigurations}

\title{
Identify distinct rate shift configurations
}
\description{
Identify topologically distinct rate shift configurations that were sampled with BAMM, and assign each sample in the posterior to one of the distinct shift configurations.	
}
\usage{
distinctShiftConfigurations(ephy, threshold)
}
\arguments{
	\item{ephy}{
		an object of class \code{bammdata}.
	}
	\item{threshold}{
		Either an object of class \code{branchprior} or a numeric value indicating the minimum frequency at which a shift has to be sampled to be included in the calculation.
	}
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
An object of class \code{branchprior} containing the marginal probability for each branch that was used to distinguish between core and non-core shifts. Alternatively, a single numeric value (we ignore all shift locations with marginal probabilities less than this value during the enumerating of distinct shift configuratins). 
}	
Results are sorted by frequency:
$frequency[1] gives the most common shift configuration sampled.
$shifts[[1]] gives the corresponding node indices for that configuration.
$samplesets[[1]] gives the indices of samples with this configuration.
}

\author{
	Dan Rabosky
}
\seealso{
	\code{\link{plot.bammshifts}}, \code{\link{credibleShiftSet}}
}
\examples{
data(whales, events.whales, prior.whales)

# Get the prior distribution on the number of shifts per branch:
bp <- getBranchShiftPriors(whales, prior.whales)


ed <- getEventData(whales, events.whales, burnin=0.25)

sc <- distinctShiftConfigurations(ed, threshold=bp)

plot.bammshifts(sc, ed, rank=1)

}



