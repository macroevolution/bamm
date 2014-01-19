\name{getmrca}
\alias{getmrca}

\title{
Find most recent common ancestors
}
\description{
Calculates the most recent common ancestor for each pair of tips. Used internally by \link{getEventData}
}
\usage{
getmrca(phy, t1, t2)
}

\arguments{
  \item{phy}{
an object of class \code{phylo}.
}
  \item{t1}{
a vector of mode integer or character corresponding to tips in \code{phy}.
}
  \item{t2}{
a vector of mode integer or character corresponding to tips in \code{phy}.
}
}
\details{
Finds the most recent common ancestor for each pair of tips where pairs are
defined as (\code{t1}[1], \code{t2}[1]), (\code{t1}[2], \code{t2}[2]), ... , (\code{t1}[i], \code{t2}[i]), ... ,(\code{t1}[n], \code{t2}[n]).
}
\value{
A vector of node numbers of the common ancestor for each pair of tips.
}
 
\author{
Mike Grundler
}
 
\keyword{ manip }
