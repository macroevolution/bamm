\name{getBranchShiftPriors}
\alias{getBranchShiftPriors}
\title{
Compute prior odds of a rate shift on each branch of a phylogeny from BAMM output
}
\description{
Computes the prior probability of a rate shift event for each branch. These results are important for identifying topological rate shift locations on phylogenies with marginal probabilities that exceed those predicted under the prior alone. 
}
\usage{
getBranchShiftPriors(phy, expectedNumberOfShifts)
}
\arguments{
  \item{phy}{
an object of class \code{phylo}
}
  \item{expectedNumberOfShifts}{
The expected number of shifts under the prior. 
}
 
}
\details{

This function computes the prior odds on the distribution of numbers of rate shift events per branch under the prior. It returns an object of class 'branchprior', which is nothing more than a copy of the original phylogenetic tree but where each branch length has been replaced by the prior probability of a rate shift on each branch.

The significance of this function is that it lets us explicitly determine which branches have shift probabilities that are elevated relative to the prior expectation.
 
}
\value{
A class \code{phylo} with all the components of the original class \code{phylo} object, with the following changes:

\item{edge.length}{
branch lengths now represent the prior probability of a rate shift on each branch
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
data(whales)
prior_tree1 <- getBranchShiftPriors(whales, expectedNumberOfShifts = 1)
prior_tree10 <- getBranchShiftPriors(whales, expectedNumberOfShifts = 10)
# plot prior expectations for branches based on these two counts:
plot(prior_tree1$edge.length ~ prior_tree10$edge.length, xlim=c(0,0.05), ylim=c(0,0.05), asp=1)
lines(x=c(0,1), y=c(0,1))
}
 
\keyword{ models }
