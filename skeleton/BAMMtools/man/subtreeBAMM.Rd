\name{subtreeBAMM}
\alias{subtreeBAMM}

\title{
Pulls out a subtree from bammdata object
}
\description{
Given a set of tips or a node, this function extracts the corresponding subtree from the \code{bammdata} object. User should specify either a set of tips or a node, and the node will overwrite the tips if both are given.
}
\usage{
subtreeBAMM(ephy,tips=NULL,node=NULL)
}

\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{tips}{
an integer or character vector indicating which tips (more than one) to be included in the subtree. 
}
  \item{node}{
an integer indicating the root of the subtree to be extracted, and it must correspond to an innernode on the tree.
}
}
\details{

This function allows users to extract a subtree from a big \code{bammdata} object, and examine the subset using \code{\link{plot.bammdata}}

}

\author{
Huateng Huang
}

\seealso{
\code{\link{getmrca}}, \code{\link{plot.bammdata}}
}
\examples{
data(whales, events.whales)
ephy <- getEventData(whales, events.whales, burnin=0.25)

#specify a set of tips for the subtree
tips <- sample(ephy$tip.label,size=20,replace=FALSE)
subphy <- subtreeBAMM(ephy,tips=tips)

#specify a innernode for subsetting
subphy <- subtreeBAMM(ephy,node=103)

#plot the subtree
plot(subphy)
}

\keyword{ graphics }