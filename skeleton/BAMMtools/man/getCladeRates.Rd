\name{getCladeRates}
\alias{getCladeRates}
\title{
Compute clade-specific mean rates
}
\description{
Computes marginal clade-specific rates of speciation, extinction, or (if relevant) trait evolution from \code{BAMM} output. 
}
\usage{
getCladeRates(ephy, node = NULL, nodetype = "include", verbose = FALSE)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
  \item{node}{
If computing rates for a specific portion of tree, the node subtending the relevant subtree
}
  \item{nodetype}{
Either "include" or "exclude". If \code{nodetype = "include"}, the rates returned by the function will be for the subtree defined by \code{node}. If \code{nodetype = "exclude"}, will compute mean rates for the tree after excluding the subtree defined by \code{node}.
}
  \item{verbose}{
If \code{TRUE}, will print the sample index as mean rates are computed for each sample from posterior. Potentially useful for extremely large trees.
}
}
\details{
Computes the time-weighted mean evolutionary rate for a given clade. Conversely, one can compute the rate for a given phylogeny while excluding a clade; this operation will give the "background" rate. It is important to understand several aspects of these mean rates. First, rates in the \code{BAMM} framework are not constant through time. Hence, the function computes the mean time-integrated rates across the subtree. Operationally, this is done by integrating the speciation rate with respect to time along each branch in the subtree. These time-integrated rates are then summed, and the sum is divided by the total sum of branch lengths for the subtree. 

The function computes a rate for each sample in the posterior, and returns a list of rate vectors. Each rate in the corresponding vector is a mean rate for a particular sample from the posterior. Hence, you can think of the return value for this function as an estimate of the marginal distribution of rates for the focal clade. You can compute means, medians, quantiles, etc from these vectors. 
}
\value{
A list with the following components:
  \item{lambda}{
A vector of speciation rates (if applicable), with the i'th rate corresponding to the mean rate from the i'th sample in the posterior
}
  \item{mu}{
A vector of extinction rates (if applicable), with the i'th rate corresponding to the mean rate from the i'th sample in the posterior
}
  \item{beta}{
A vector of phenotypic rates (if applicable), with the i'th rate corresponding to the mean rate from the i'th sample in the posterior
}
}
\references{
\code{bamm-project.org}
}
\author{
Dan Rabosky
}
 
\examples{
data(events.whales)
data(whales)
ed <- getEventData(whales, events.whales)
all_rates <- getCladeRates(ed)

mean(all_rates$lambda)
mean(all_rates$mu)
# joint density of mean speciation and extinction rates:
plot(all_rates$mu ~ all_rates$lambda)

# clade specific rates: here for Dolphin subtree:
dol_rates <- getCladeRates(ed, node=140)
mean(dol_rates$lambda)
mean(dol_rates$mu)

}
 
\keyword{ models }
