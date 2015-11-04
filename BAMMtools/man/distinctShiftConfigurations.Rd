\name{distinctShiftConfigurations}
\alias{distinctShiftConfigurations}

\title{
Identify distinct rate shift configurations
}
\description{
Identify topologically distinct rate shift configurations that were sampled with BAMM, and assign each sample in the posterior to one of the distinct shift configurations.	
}
\usage{
distinctShiftConfigurations(ephy, expectedNumberOfShifts, threshold, ...)
}
\arguments{
	\item{ephy}{
		an object of class \code{bammdata}.
	}
	\item{expectedNumberOfShifts}{
		The expected number of rate shifts under the prior
	}
	
	\item{threshold}{
		Threshold value for marginal posterior-to-prior odds ratios, used to identify branches with elevated shift probabilities relative to prior (core vs non-core shifts)
	}
	\item{...}{
		Other arguments to distinctShiftConfigurations (possibly deprecated args)
	}
}
\details{
See Rabosky et al (2014) and the BAMM project website for details on the nature of distinct shift configurations and especially the distinction between "core" and "non-core" rate shifts. Note that branches with elevated marginal posterior probabilities relative to the prior (marginal odds ratios) cannot be claimed to have "significant" evidence for a rate shift on the basis of this evidence alone.

}


\value{
an object of class \code{bammshifts}. This is a list with the following components
\item{marg.probs}{
A list of the marginal probability of a shift occurring at each node of the phylogeny for each
distinct rate shift configuration.	
}
\item{marginal_odd_ratio}{
Marginal posterior-to-prior odds ratios for one or more rate shifts an a given branch
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
\item{coreshifts}{
A vector of node numbers corresponding to the core shifts. All of these nodes have a marginal odds ratio of at least \code{threshold} supporting a rate shift.
}
\item{threshold}{
A single numeric value giving the marginal posterior:prior odds ratio threshold used during enumeration of distinct shift configurations
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
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)

sc <- distinctShiftConfigurations(ed, expectedNumberOfShifts = 1, threshold = 5)

plot.bammshifts(sc, ed, rank=1)

}



