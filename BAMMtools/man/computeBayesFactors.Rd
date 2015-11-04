\name{computeBayesFactors}
\alias{computeBayesFactors}
\title{
computeBayesFactors
}
\description{
Computes pairwise Bayes factors for a set of macroevolutionary models sampled using \code{BAMM}, using MCMC simulation output.
}
\usage{
computeBayesFactors(postdata, expectedNumberOfShifts, burnin = 0.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{postdata}{
Filename for the MCMC output file from a \code{BAMM} run. Alternatively, a dataframe containing this information.
}
  \item{expectedNumberOfShifts}{
Expected number of shifts under the prior
}
  \item{burnin}{
What fraction of samples to discard from postdata as burnin?
}
 \item{...}{
 Additional arguments to computeBayesFactors
 }
}
\details{
This function returns a matrix of pairwise Bayes factors, where the Bayes factor is the ratio of marginal likelihoods between two models M_{i} and M_{j}. Numerator models are given as rows, and denominator models as columns. Row names and column names give the number of shifts in the corresponding model. Suppose you have an output matrix with row and column names 0:3 (0, 1, 2, 3). Model 0 is a model with just a single process (starting at the root), and no among-lineage rate heterogeneity. 

If \code{computeBayesFactors} gives a matrix \code{mm}, and \code{mm[2,1]} is 10.0, this implies Bayes factor evidence of 10 in favor of the 2nd row model (a model with 1 process; e.g., \code{rownames(mm)[2]} over the first column model (a model with a single process).

This function will only compute Bayes factors between models which were actually sampled during simulation of the posterior. Hence, if a model has such low probability that it is never visited by BAMM during the simulation of the posterior, it will be impossible to estimate its posterior probability (and thus, you will get no Bayes factors involving this particular model). This is likely to change in the future with more robust methods for estimating posterior probabilities in the tails of the distribution.  
}
\value{
A matrix of pairwise Bayes factors between models.
}

\author{

Dan Rabosky

}
\examples{
data(mcmc.whales)
 
computeBayesFactors(mcmc.whales, expectedNumberOfShifts = 1, burnin=0.1)
}

\keyword{ models }
