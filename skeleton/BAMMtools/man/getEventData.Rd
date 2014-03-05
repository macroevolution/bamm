\name{getEventData}
\alias{getEventData}
\title{
Create \code{bammdata} object from MCMC output
}
\description{
\code{getEventData} reads shift configuration data (the "event data" output) from a \code{BAMM} analysis and creates a \code{bammdata} object. The \code{bammdata} object is fundamental for extracting information about macroevolutionary rate variation through time and among lineages.
}
\usage{
getEventData(phy, eventdata, burnin = 0, nsamples = NULL, 
	verbose = FALSE, type = "diversification")
}
\arguments{
  \item{phy}{
an object of class \code{phylo} - specifically, the time calibrated tree that was analyzed with \code{BAMM}.
}
  \item{eventdata}{
a character string specifying the path to a \code{BAMM} event-data file. Alternatively, an object of class \code{data.frame} that includes the event data from a \code{BAMM} run.
}
  \item{burnin}{
a numeric indicating the fraction of posterior samples to discard as burn-in.
}
  \item{nsamples}{
an integer indicating the number of posterior samples to include in the \code{bammdata} object. May be \code{NULL}.
}
  \item{verbose}{
a logical. If \code{TRUE} progess is outputted to the console. Defaults to \code{FALSE}.
}
  \item{type}{
a character string. Either "diversification" or "trait" depending on your BAMM analysis.
}
}

\details{
In the BAMM framework, an "event" defines a macroevolutionary process of diversification 
or trait evolution. Every sample from the posterior includes at least one process, defined by such an "event". If a given sample includes just a single event, then the dynamics of diversification or trait evolution can be described entirely by a single time-constant or time-varying process that begins at the root of the tree. Any sample from the posterior distribution may include a complex mixture of distinct processes. To represent temporal heterogeneity in macroevolutionary rates, \code{BAMM} models a rate \eqn{R}, e.g. speciation, as a function that changes exponentially with time: \eqn{R(t) = R(0)*exp(b*t)}. Here \eqn{R(0)} is the initial rate and \eqn{b} is a parameter determining how quickly that rate grows or decays with time. 

The \code{eventdata} file (or data frame) is a record of events and associated parameters that were sampled with \code{BAMM} during simulation of the posterior with reversible jump MCMC. This complex, information-rich file is processed into a \code{bammdata} object, which serves as the core data object for numerous downstream analyses. From a \code{bammdata} object, you can summarize rate variation through time, among clades, extract locations of rate shifts, summarize clade-specific rates of speciation and extinction, and more.

In general, the user does not need to be concerned with the details of a \code{bammdata} object. The object is used as input by a number of \code{BAMMtools} functions. 
 
The parameter \code{nsamples} can be used to reduce the total amount of data included in the raw eventdata output from a \code{BAMM} run. The final \code{bammdata} object will consist of all data for \code{nsamples} from the posterior. These \code{nsamples} are equally spaced after discarding some \code{burnin} fraction as "burn-in". If \code{nsamples} is set to \code{NULL}, the \code{bammdata} object will include all samples in the posterior after discarding the \code{burnin} fraction.

}
\value{
A list with many components:
\item{edge}{See documentation for class \code{phylo} in package ape.}
\item{Nnode}{See documentation for class \code{phylo} in package ape.}
\item{tip.label}{See documentation for class \code{phylo} in package ape.}
\item{edge.length}{See documentation for class \code{phylo} in package ape.}

\item{begin}{The beginning time of each branch in absolute time (the root is set to time zero)}
\item{end}{The ending time of each branch in absolute time.}

\item{numberEvents}{
An integer vector with the number of events contained in \code{phy} for each posterior sample. The length of this vector is equal to the number of posterior samples in the \code{bammdata} object.
}
\item{eventData}{
A list of dataframes. Each element is a single posterior sample.
Each row in a dataframe holds the data for a single event.
Data associated with an event are: \code{node} - a node number. This identifies 
the branch where the event originates. \code{time} - this is the absolute time on 
that branch where the event originates (with the root at time 0). \code{lam1} - an initial rate of speciation or trait evolution.
\code{lam2} - a decay/growth parameter. \code{mu1} - an initial rate of extinction.
\code{mu2} - a decay/growth parameter. \code{index} - a unique integer associated with the event.
See 'Details'.
}
\item{eventVectors}{
A list of integer vectors. Each element is a single posterior sample.
For each branch in \code{phy} the index of the event that occurs along that branch.
Branches are ordered increasing here and elsewhere.
}
\item{eventBranchSegs}{
A list of matrices. Each element is a single posterior sample.
Each matrix has four columns: \code{Column 1} identifies a node in \code{phy}. 
\code{Column 2} identifies the beginning time of the branch or segment 
of the branch that subtends the node in \code{Column 1}. \code{Column 3} identifies the 
ending time of the branch or segment of the branch that subtends the node in \code{Column 1}. 
\code{Column 4} identifies the index of the event that occurs along the branch or segment 
of the branch that subtends the node in \code{Column 1}.
}
\item{tipStates}{
A list of integer vectors. Each element is a single posterior sample.
For each tip the index of the event that occurs along the branch subtending the tip.
Tips are ordered increasing here and elsewhere.
}
\item{tipLambda}{
A list of numeric vectors. Each element is a single posterior sample.
For each tip the rate of speciation or trait evolution at the end of 
the terminal branch subtending that tip.
}
\item{tipMu}{
A list of numeric vectors. Each element is a single posterior sample.
For each tip the rate of extinction at the end of the terminal branch
subtending that tip. Meaningless if working with \code{BAMM} trait results.
}
\item{meanTipLambda}{
For each tip the mean of the marginal posterior density of the rate of speciation
or trait evolution at the end of the terminal branch subtending that tip.
}
\item{meanTipMu}{
For each tip the mean of the marginal posterior density of the rate of extinction
at the end of the terminal branch subtending that tip. Meaningless if working with BAMM trait results.
}
\item{type}{
A character string. Either "diversification" or "trait" depending on your BAMM analysis.
}
\item{downseq}{
An integer vector holding the nodes of \code{phy}. The order corresponds to the 
order in which nodes are visited by a pre-order tree traversal.
}
\item{lastvisit}{
An integer vector giving the index of the last node visited by the node in the corresponding
position in \code{downseq}. \code{downseq} and \code{lastvisit} can be used to quickly retrieve
the descendants of any node. e.g. the descendants of node 89 can be found by 
\code{downseq[which(downseq==89):which(downseq==lastvisit[89])}.
}

}
\references{
\url{bamm-project.org}

Rabosky, D. L. 2014. Automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees. PLoS ONE 9:e89543

Rabosky, D. L., F. Santini, J. T. Eastman, S. A. Smith, B. L. Sidlauskas, J. Chang, and M. E. Alfaro. 2013. Rates of speciation and morphological evolution are correlated across the largest vertebrate radiation. Nature Communications DOI: 10.1038/ncomms2958 

}
\author{
Dan Rabosky, Mike Grundler
}
\note{
Currently the function does not check for duplicate tip labels in \code{phy}, which may cause the function to choke.
}

\seealso{
\code{\link{summary.bammdata}}, \code{\link{plot.bammdata}}, \code{\link{dtRates}}.
}
\examples{
data(primates)
data(events.primates)
xx <- getEventData(primates, events.primates, burnin=0.25, type = 'trait')

# compute mean phenotypic rate for primate body size evolution:
brates <- getCladeRates(xx)
mean(brates$beta)

# Plot rates:
plot.bammdata(xx)

}

\keyword{ graphics }
\keyword{ models }

