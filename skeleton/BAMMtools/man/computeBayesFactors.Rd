\name{computeBayesFactors}
\alias{computeBayesFactors}
\title{
computeBayesFactors
}
\description{
Computes pairwise Bayes factors for a set of macroevolutionary models sampled using \code{BAMM}, using MCMC simulation output.
}
\usage{
computeBayesFactors(postdata, priordata, burnin = 0.1, modelset = NULL, 
	threshpost = 1, threshprior = 0, nbprior = FALSE, strict = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{postdata}{
Filename for the MCMC output file from a \code{BAMM} run. Alternatively, a dataframe containing this information.
}
  \item{priordata}{
Filename for the \code{BAMM}-generated output file containing a simulation of shifts under the prior. Alternatively, a dataframe containing this information.
}
  \item{burnin}{
What fraction of samples to discard from postdata as burnin?
}
  \item{modelset}{
Integer set of models for which you wish to compute Bayes factors. An integer value of k implies a model with k rate shifts.
}
  \item{threshpost}{
Only compute Bayes factors for models sampled at least this many times in the posterior.
}
  \item{threshprior}{
Only compute Bayes factors for models sampled at least this many times during simulation of the prior.
}
\item{nbprior}{
Use a negative binomial distribution to approximate the prior.	
}
  \item{strict}{
If TRUE, requires that both models i and j be sampled at least \code{threshpost} and \code{threshprior} times. Otherwise, only requires that one of the target models meet this condition to be included in the calculation.
}
}
\details{
This function returns a matrix of pairwise Bayes factors, where the Bayes factor is the ratio of marginal likelihoods between two models M_{i} and M_{j}. Numerator models are given as rows, and denominator models as columns. Row names and column names give the number of shifts in the corresponding model. Suppose you have an output matrix with row and column names 0:3 (0, 1, 2, 3). Model 0 is a model with just a single process (starting at the root), and no among-lineage rate heterogeneity. 

If \code{computeBayesFactors} gives a matrix \code{mm}, and \code{mm[2,1]} is 10.0, this implies Bayes factor evidence of 10 in favor of the 2nd row model (a model with 1 process; e.g., \code{rownames(mm)[2]} over the first column model (a model with a single process).


}
\value{
A matrix of pairwise Bayes factors between models.
}

\author{

Dan Rabosky

}

\examples{
data(mcmc.whales)
data(prior.whales)
computeBayesFactors(mcmc.whales, prior.whales, burnin=0.1)
}

\keyword{ models }
