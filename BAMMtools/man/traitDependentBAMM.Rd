\name{traitDependentBAMM}
\alias{traitDependentBAMM}

\title{
STRAPP: STructured Rate Permutations on Phylogenies
}
\description{
Given a \code{bammdata} object and a vector of trait data (continuous), assess whether the correlation between the trait and bamm estimated speciation , extinction or net diversification rate is significant using permutation. A set of posterior samples is randomly drawn from the \code{bammdata} object. If the trait is continuous, this function calculates the correlation coefficients between the trait and tip rates (observed correlation), as well as that with permuted rates for each posterior sample. In a one-tailed test for positive correlations, the reported p-value is the proportion of the posterior samples in which the observed correlation is larger than the correlations calculated with permuted rates. In a two-tailed test, the reported p-value is the proportion of the posterior samples in which the null correlation is as extreme as the correlation observed. If the trait is binary, the U statistic of the Mann-Whitney test is calculated instead of correlation coefficients to assess whether there is a significant difference in rate between the two trait states. For categorical traits with more than two states, the Kruskal-Wallis rank sum statistic is used.
}
\usage{
traitDependentBAMM(ephy, traits, reps, rate = 'speciation', return.full = FALSE, 
  method = 'spearman', logrates = TRUE, two.tailed = TRUE, 
  traitorder = NA, nthreads = 1) 
}

\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{traits}{
a vector of trait data, with names corresponding to tips in the \code{bammdata} object. It can be numeric or categorical. 
}
  \item{reps}{
an integer specifying the number of permutations (i.e., the number of posterior samples to randomly draw with replacement from the \code{bammdata} object)
}
  \item{rate}{
a character string specifying which estimated rate from the \code{bammdata} object to use for testing correlation, must be one of 'speciation', 'extinction', or 'net diversification'. Defaults to 'speciation'. You can specify just the initial letter.
}
  \item{return.full}{
a logical. If \code{TRUE}, the list of posterior samples, the observed correlation for each posterior sample, and the null distribution will be included in the returned object. Defaults to \code{FALSE}.
}
  \item{method}{
a character string, must be one of 'spearman', 'pearson', "mann-whitney", or "kruskal". Defaults to 'spearman'. You can specify just the initial letter.
}
  \item{logrates}{
a logical. If \code{TRUE} log-transform the rates before analysis. Defaults to \code{TRUE}. This can only matter for the pearson correlation. 
}
  \item{two.tailed}{
a logical, used for cotinuous trait data. If \code{TRUE}, perform a two tailed statistical test (i.e., if the null distribution is symmetric, it is equivalent to doubling the p-value). Defaults to \code{TRUE}.  
}
  \item{traitorder}{
a character string specifying the direction of correlation for the alternative hypothesis. For continuous traits, it must be either "positive" or "negative"; only the initial letter is needed. For binary traits, it must be a string indicating states with increasing rate under the alternative hypothesis, separated by comma (e.g., 'A, B'). One-tailed test for categorical data with more than two states is not supported.
  }
  	\item{nthreads}{
number of threads to use for parallelization of the function. The R package \code{snow} must be loaded for \code{nthreads > 1}.
  	}
}

\details{
Tip rates-- speciation, extinction or net diversification rates-- are permuted in a way such that pairwise covariances in rates between species are maintained. That is, tips with the same \code{tipStates} still have the same rate after permutation. Posterior samples are randomly selected with replacement from the \code{bammdata} object, so \code{reps} could be smaller or larger than the total number of samples in the object.
}
\value{
  \item{estimate}{
a numeric value for continous trait data- the average observed correlation between tip rates and the trait across the posterior samples. For categorical traits, it is a list showing the median species-specific rates for each trait state. 
}
  \item{p.value}{
a numeric value. The probability that the observed correlation is less than or equal to a sample from the null distribution.
}
  \item{method}{
a character string, as input.
}
  \item{rate}{
a character string, as input.
}
  \item{two.tailed}{
a logical, as input.
}
  \item{gen}{
an integer vector, recording which posterior samples were selected. Only present when \code{return.full} is \code{TRUE}.
}
  \item{obs.corr}{
a numeric vector, the observed correlation coefficents for each posterior sample. Only present when \code{return.full} is \code{TRUE}. For binary traits, centered U statistics (U - n1* n2/2; where n1 and n2 are the number of species in each state of the binary trait) is reported
}
  \item{null}{
a numeric vector. The null distribution of correlation coefficients (or centered U statistics for binary traits) from permutation. Only present when \code{return.full} is \code{TRUE}.
}
}
\references{
\url{bamm-project.org}

Rabosky, D. L. 2014. Automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees. PLoS ONE 9:e89543

Rabosky, D. L., F. Santini, J. T. Eastman, S. A. Smith, B. L. Sidlauskas, J. Chang, and M. E. Alfaro. 2013. Rates of speciation and morphological evolution are correlated across the largest vertebrate radiation. Nature Communications DOI: 10.1038/ncomms2958 

}
\author{
Daniel L. Rabosky, Huateng Huang
}

\seealso{
\code{\link{subtreeBAMM}}
}
\examples{
# using a small subset of the fish data set (300 taxa) in Rabosky et al. 2013. Nat. Com. paper
data(fishes)
data(events.fishes)
xx <- getEventData(phy = fishes, eventdata = events.fishes, type = "diversification")
#traits.fishes is the trait -- body size
data(traits.fishes)
x <- traitDependentBAMM(ephy = xx, traits = traits.fishes, reps = 1000, return.full = TRUE, 
		method = 's', logrates = TRUE, two.tailed = TRUE)

}
\keyword{ nonparametric }
\keyword{ univar }% __ONLY ONE__ keyword per line
