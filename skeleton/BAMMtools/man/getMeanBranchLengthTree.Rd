\name{getMeanBranchLengthTree}
\alias{getMeanBranchLengthTree}
\title{
Compute phylogeny with branch lengths equal to corresponding macroevolutionary rate estimates
}
\description{
Takes a \code{bammdata} object and computes a phylogenetic tree where branch lengths are equal to the mean of the marginal distributions of rates on each branch. This tree can be plotted to visualize rate variation across a phylogeny. 
}
\usage{
getMeanBranchLengthTree(ephy, rate = "speciation")
}

\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{rate}{
the type of rate-tree to be computed. Options: "speciation" (default), "extinction", "ndr" (net diversification), and "trait"
}
}
\value{
A list with the following components:
  \item{phy}{
A phylogenetic tree, topologically identical to the model tree, but with branch lengths replaced by the mean (marginal) rates on each branch as estimated from the posterior samples in the \code{bammdata} object.
}
  \item{mean}{
The mean rate over all branches
}
  \item{median}{
the median rate over all branches
}

}
\references{
\url{www.bamm-project.org} 
}
\author{
Dan Rabosky
}

\seealso{
\link{plot.bammdata}
}
\examples{
data(whales)
data(events.whales)
ed <- getEventData(whales, events.whales)
ed2 <- subsetEventData(ed, index = 1:20)
ratetree <- getMeanBranchLengthTree(ed2, rate='speciation')
plot(ratetree$phy, show.tip.label=FALSE)
}
\keyword{ models}
