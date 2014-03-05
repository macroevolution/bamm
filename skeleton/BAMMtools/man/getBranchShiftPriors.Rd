\name{getBranchShiftPriors}
\alias{getBranchShiftPriors}
\title{
Compute expected percentiles of distribution of rate shifts per branch under prior
}
\description{
Computes the 95th percentile of the distribution of rate shift events per branch under the prior. These results are important for identifying topological rate shift locations on phylogenies with marginal probabilities that exceed those predicted under the prior alone. 
}
\usage{
getBranchShiftPriors(phy, priordata)
}
\arguments{
  \item{phy}{
an object of class \code{phylo}
}
  \item{priordata}{
Filename for the \code{BAMM}-generated output file containing a simulation of shifts under the prior. Alternatively, a dataframe containing this information. 
}
}
\details{

This function computes the 95th percentile on the distribution of numbers of rate shift events per branch under the prior. It returns an object of class 'branchprior', which is nothing more than a copy of the original phylogenetic tree but where each branch length has been replaced by the number of 95th percentile on the distribution of shifts per branch. 

Thus, a branch with a value of 0.05 means that, on average, we expect (with probability 0.95) no more than 0.05 shift events on this branch under the prior alone. 

The significance of this function is that it lets us explicitly determine which branches have shift probabilities that are elevated relative to the prior expectation. The user does not work with the output of this function directly; it is passed to other functions \link[BAMMtools]{credibleShiftSet} and \link[BAMMtools]{distinctShiftConfigurations}.
 
}
\value{
A class \code{branchprior} with all the components of the original class \code{phylo} object, with the following changes:

\item{edge.length}{
branch lengths now represent the expected number of shifts per branch at the 95th percentile under the prior. 
}
 
\item{criterion}{
the relevant percentile (usually 0.95)
}
}
\references{
\code{bamm-project.org}
}
\author{
Dan Rabosky
}

\seealso{
\code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammshifts}}, \code{\link{summary.credibleshiftset}}, \code{\link{plot.credibleshiftset}}, \code{\link{credibleShiftSet}}	
} 
 
\examples{
data(events.whales)
data(whales)
data(prior.whales)
ed <- getEventData(whales, events.whales)

# Get branch-specific cutoffs for identifying branches with 
#	elevated shift probabilities:
bprior <- getBranchShiftPriors(whales, prior.whales)
cset <- credibleShiftSet(ed, threshold = bprior)

# Here is the total number of samples in the posterior:
length(ed$eventData)

# And here is the number of distinct shift configurations:
cset$number.distinct

# here is the summary statistics:
summary(cset)

# Accessing the raw frequency vector for the credible set:
cset$frequency

#The cumulative frequencies:
cset$cumulative

# The first element is the shift configuration with the maximum 
#  a posteriori probability. We can identify all the samples from 
#  posterior that show this shift configuration:

cset$indices[[1]]
 
 
 
}
 
\keyword{ models }
