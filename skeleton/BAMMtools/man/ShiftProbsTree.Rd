\name{ShiftProbsTree}
\alias{marginalShiftProbsTree}
\alias{cumulativeShiftProbsTree}

\title{
Branch-specific rate shift probabilities
}

\description{
\code{marginalShiftProbsTree} computes a version of a phylogenetic tree where each branch length is equal to the marginal probability that a shift occurred on a particular branch. The \code{cumulativeShiftProbsTree} includes the cumulative probability that a shift occurred on a given branch. See details.

}

\usage{
cumulativeShiftProbsTree(ephy)
marginalShiftProbsTree(ephy)
}

\arguments{
  \item{ephy}{
An object of class \code{bammdata}
}

}
\details{

The \emph{marginal shift probability tree} is a copy of the target phylogeny, but where each branch length is equal to the branch-specific marginal probability that a rate-shift occurred on the focal branch. For example, a branch length of 0.333 implies that 1/3 of all samples from the posterior had a rate shift on the focal branch. 

\bold{Note:} It is highly inaccurate to use marginal shift probabilities as a measure of whether diversification rate heterogeneity occurs within a given dataset. Consider the following example. Suppose you have a tree with topology (A, (B, C)). You find a marginal shift probability of 0.5 on the branch leading to clade C, and also a marginal shift probability of 0.5 on the branch leading to clade BC. Even though the marginal shift probabilities appear low, it may be the case that the joint probability of a shift occurring on \emph{either} the branch leading to C or BC is 1.0. Hence, you could be extremely confident (posterior probabilities approaching 1.0) in rate heterogeneity, yet find that no single branch has a particularly high marginal shift probability. In fact, this is exactly what we expect in most real datasets, because there is rarely enough signal to strongly support the occurrence of a shift on any particular branch.

The \emph{cumulative shift probability tree} is a copy of the target phylogeny but where branch lengths are equal to the cumulative probability that a rate shift occurred somewhere on the path between the root and the focal branch. A branch length equal to 0.0 implies that the branch in question has evolutionary rate dynamics that are shared with the evolutionary process starting at the root of the tree. A branch length of 1.0 implies that, with posterior probability 1.0, the rate dynamics on a branch are decoupled from the "root process". 

}
\value{
An object of class \code{phylo}, but with branch lengths equal to the marginal or cumulative shift probabilities.
}
\references{
\url{bamm-project.org}
}
\author{
Dan Rabosky
}

\seealso{
\link[BAMMtools]{maximumShiftCredibility}
}
\examples{
data(whales)
data(events.whales)
ed <- getEventData(whales, events.whales)

# computing the marginal shift probs tree:
mst <- marginalShiftProbsTree(ed)

# The cumulative shift probs tree:
cst <- cumulativeShiftProbsTree(ed)

#compare the two types of shift trees side-by-side:
plot.new()
par(mfrow=c(1,2))
plot.phylo(mst, no.margin=TRUE, show.tip.label=FALSE)
plot.phylo(cst, no.margin=TRUE, show.tip.label=FALSE)
}
 
\keyword{ graphics}
