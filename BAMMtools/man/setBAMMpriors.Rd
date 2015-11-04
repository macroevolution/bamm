\name{setBAMMpriors}
\alias{setBAMMpriors}

\title{
setBAMMpriors
}
\description{
Set priors for BAMM analysis

}
\usage{
setBAMMpriors(phy, total.taxa = NULL, traits=NULL, outfile = 'myPriors.txt', 
	Nmax = 1000, suppressWarning = FALSE)
}

\arguments{
  \item{phy}{
an object of class \code{phylo}, e.g., the phylogenetic tree you will analyze with BAMM
}
  \item{total.taxa}{
If doing speciation-extinction analysis, the total number of taxa in your clade. If your tree contains all taxa in the clade (100\% sampling), then leave this as NULL
}
  \item{traits}{
A filename where the trait data (BAMM format) are stored. Leave NULL if doing a speciation-extinction analysis.
}
  \item{outfile}{
Filename for outputting the sample priors block. If \code{NULL}, then a vector is returned instead.
}
  \item{Nmax}{
If analyzing a very large tree for phenotypic evolution, uses only this many taxa to estimate priors for your dataset. Avoid matrix inversion issues with large numbers of tips. 
}
	\item{suppressWarning}{
if \code{TRUE}, then the warning about setting the Poisson rate prior is suppressed. Only applies if \code{outfile = NULL}.
	}
}

\details{
This is a "quick and dirty" tool for identifying approximately acceptable priors for a BAMM analysis. We have found that choice of prior can have a substantial impact on BAMM analyses. It is difficult to simply set a default prior that applies across datasets, because users often have trees with branch lengths in very different units (e.g., numbers of substitutions versus millions of years). Hence, without some careful attention, you can inadvertently specify some very bad prior distributions. This function is designed to at least put you in the right ballpark for decent prior distributions, but there are no guarantees that these are most appropriate for your data.

The general rules applied here are: 

For the \code{lambdaInitPrior}, we estimate the speciation rate of the data under a pure birth model. We then set this prior to give an exponential distribution with a mean five times greater than this computed pure birth speciation estimate.
 
The \code{lambdaShiftPrior} is the standard deviation of the normal prior on the exponential change parameter k. We set the prior distribution based on the age of the root of the tree. We set the standard deviation of this distribution such that 2 standard deviations give a parameter that will yield a 90\% decline in the initial speciation rate between the root of the tree and the tips. The basic model is lambda(t) = lambda_0 * exp(k * t). This is a straightforward calculation: let x = -log(0.1) / TMAX, where TMAX is the age of the tree. Then set the standard deviation equal to (x / 2).

We set \code{muInitPrior} equal to \code{lambdaInitPrior}.  

For trait evolution, we first compute the maximum likelihood estimate of the variance of a Brownian motion process of trait evolution. The prior \code{betaInitPrior} is then set to an exponential distribution with a mean 5 times greater than this value (similar to what is done for lambda and mu, above).  

This function generates an output file containing a prior parameters block that can be pasted directly into the priors section of your \code{BAMM} control file.

}

\value{
The function does not return anything. It simply performs some calculations and writes formatted output to a file.
However, if \code{outfile = NULL}, then a named vector is returned. 

}
 
\author{
Dan Rabosky
}

\examples{
data(whales)

setBAMMpriors(phy = whales, total.taxa = 89, outfile = NULL)


}
 
 
\keyword{ models}
