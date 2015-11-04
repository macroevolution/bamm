\name{getMarginalBranchRateMatrix}
\alias{getMarginalBranchRateMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute mean branch rates for \code{bammdata} object 
}
\description{
For each sample in the posterior, computes the mean rate for each branch in the focal phylogeny (speciation, extinction, trait evolution). If the \code{bammdata} object contains \emph{nsamples} samples and the target phylogeny has \emph{nbranches} branches, the function will compute a matrix of \emph{nbranches} x \emph{nsamples}.
}
\usage{
getMarginalBranchRateMatrix(ephy, verbose = FALSE)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{verbose}{
print progress during processing of \code{bammdata} object
}
}
\details{
If a \code{type = 'diversification'} \code{bammdata} object is passed as an argument, the function will return matrices for both speciation and extinction. If \code{type = 'trait'} object, the matrix will simply be the corresponding phenotypic rates. Branch-specific rates are the mean rates computed by integrating the relevant rate-through-time function along each branch, then dividing by the length of the branch.

}
\value{
Returns a list with the following components:
\item{lambda_branch_matrix}{
A \code{nbranches x nsamples} matrix of mean speciation rates for each branch.
}
\item{mu_branch_matrix}{
A \code{nbranches x nsamples} matrix of mean extinction rates for each branch.

}
\item{beta_branch_matrix}{
A \code{nbranches x nsamples} matrix of mean phenotypic rates for each branch.

}

}
 
\author{
Dan Rabosky
}
 
 
\examples{
data(whales)
data(events.whales)
ed <- getEventData(whales, events.whales, nsamples = 10)
mbr <- getMarginalBranchRateMatrix(ed)
dim(mbr$lambda_branch_matrix)

}

\keyword{ models }
