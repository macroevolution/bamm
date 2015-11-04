\name{marginalOddsRatioBranches}
\alias{marginalOddsRatioBranches}
\title{
Ratio of (marginal) posterior-to-prior probabilities on individual branches
}
\description{
Compute marginal posterior-to-prior odds ratio associated with observing one or more rate shift on a given branch.
}
\usage{
marginalOddsRatioBranches(ephy, expectedNumberOfShifts)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
 \item{expectedNumberOfShifts}{
Expected number of shifts under the prior alone
 }
}
\details{
This function returns a copy of a phylogenetic tree where each branch length is equal to the marginal odds ratio in favor of a rate shift on a particular branch. These cannot be interpreted as evidence for a rate shift in an absolute sense. As explained on the website, they are a marginal odds ratio. This function is provided primarily for the purpose of distinguishing core and non-core shifts.
 
}

\value{
A object of class \code{phylo} but where each branch length is equal to the marginal shift  odds on each branch.
}

\author{
Dan Rabosky
}

\seealso{
\code{\link{getBranchShiftPriors}}, \code{\link{distinctShiftConfigurations}}, \code{\link{credibleShiftSet}}
}
\examples{
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin = 0.1)
bftree <- marginalOddsRatioBranches(ed, expectedNumberOfShifts = 1)
margprobs <- marginalShiftProbsTree(ed)
plot.new()
par(mfcol = c(1,3))
# Original tree:
plot.phylo(whales, show.tip.label = FALSE)
mtext(side = 3, text = "original tree", cex = 1)

plot.phylo(margprobs, show.tip.label = FALSE)
mtext(side = 3, text = "Marginal probs", cex = 1)

plot.phylo(bftree, show.tip.label = FALSE)
mtext(side = 3, text = "Marginal odds ratio", cex = 1)


}
