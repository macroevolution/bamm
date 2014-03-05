\name{getCohortMatrix}
\alias{getCohortMatrix}
\alias{getBAMMCorrelationMatrix}
\alias{getMacroevolutionaryCohortMatrix}
\title{
Compute the pairwise correlation in rate regimes between all tips in a bammdata object
}
\description{
Takes a \code{bammdata} object and computes the pairwise correlation in evolutionary rate regimes between all tips in the phylogeny. This can be used to identify cohorts of taxa that share common macroevolutionary rate parameters. It can also be used to construct a correlation matrix for GLS analyses using BAMM-estimated tip rates of speciation, extinction, or phenotypic evolution.
 
}
\usage{
getCohortMatrix(ephy)
}

\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
}

\details{
The cohort matrix is important for interpreting and visualizing macroevolutionary dynamics. Each entry [i, j] of the cohort matrix is the probability that taxon i and taxon j share a common macroevolutionary rate regime. To compute this, we simply tabulate the percentage of samples from the posterior where taxon i and taxon j were placed in the same rate regime. If there is no rate heterogeneity in the dataset (e.g., the data are best explained by a single rate regime), then all species will tend to share the same rate regime and all values of the cohort matrix will approach 1. 

A value of 0 between any two taxa means that at least one rate shift occurred on the nodal path connecting them in 100\% of samples from the posterior. A value of 0.50 would imply that 50\% of samples from the posterior included a rate shift on the path connecting taxa i and j. See below (examples) for an illustration of this. 

}
\value{
A numeric matrix of dimension k x k, where k is the number of species in the phylogeny included in the \code{bammdata} object. Species names are included as row names and column names. The matrix is symmetric, such that the values for entry [i , j] will mirror those for [j , i].

}

\references{
\url{www.bamm-project.org} 
}
\author{
Dan Rabosky
}


\examples{
data(whales)
data(events.whales)
ed <- getEventData(whales, events.whales)

cormat <- getCohortMatrix(ed)

dim(cormat)
hist(cormat, breaks=50)

}
\keyword{ models}
