\name{getEventData}
\alias{getEventData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load BAMM data
}
\description{
\code{getEventData} reads BAMM data and creates a \code{bammdata} object.
}
\usage{
getEventData(phy, eventdata, burnin = 0, nsamples = NULL, verbose = FALSE, type = "diversification")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
an object of class \code{phylo}.
}
  \item{eventdata}{
a character string with the path to a BAMM event file or an object of class \code{dataframe} with the 
same information.
}
  \item{burnin}{
a numeric indicating the fraction of posterior samples to discard as burn-in.
}
  \item{nsamples}{
an integer indicating the number of posterior samples to load. May be \code{NULL}.
}
  \item{verbose}{
a logical. If \code{TRUE} progess is outputted to the console. Defaults to \code{FALSE}.
}
  \item{type}{
a character string. Either "diversification" or "trait" depending on your BAMM analysis.
}
}

\details{
In the BAMM framework, an "event" is a macroevolutionary process of diversification 
or trait evolution. Every phylogenetic tree contains at least 1 event that begins 
at the root node. A macroevolutionary process of diversification comprises rates of 
speciation and extinction that vary through time and among lineages. A macroevolutionary 
process of trait evolution comprises rates of phenotypic change that vary through time 
and among lineages. 

To represent temporally heterogeneity in macroevolutionary rates
BAMM uses the exponential change function ...

To represent among lineage heterogeneity BAMM allows different parts of the tree to
evolve under independent macroevolutionary processes. The number and configuration of
processes are estimated using reversible jump MCMC ...

 ... further details ...
}
\value{
A list with many components:
\item{edge}{See documentation for class \code{phylo} in package ape.}
\item{Nnode}{See documentation for class \code{phylo} in package ape.}
\item{tip.label}{See documentation for class \code{phylo} in package ape.}
\item{edge.length}{See documentation for class \code{phylo} in package ape.}

\item{begin}{The beginning time of each branch. The root is set to time zero.}
\item{end}{The ending time of each branch.}

\item{numberEvents}{
An integer vector with the number of events contained in \code{phy} for each posterior sample.
}
\item{eventData}{
A list of dataframes. Each element is a single posterior sample.
Each row in a dataframe holds the data for a single event.
Data associated with an event are: \code{node} - a node number. This identifies 
the branch where the event originates. \code{time} - this is the absolute time on 
that branch where the event originates. \code{lam1} - an initial rate of speciation or trait evolution.
\code{lam2} - a decay/growth parameter. \code{mu1} - an initial rate of extinction.
\code{mu2} - a decay/growth parameter. \code{index} - a unique integer associated with the event.
See 'Details'.
}
\item{eventVectors}{
A list of integer vectors. Each element is a single posterior sample.
For each branch in \code{phy} the index of the event that obtains along that branch.
Branches are ordered increasing here and elsewhere.
}
\item{eventBranchSegs}{
A list of matrices. Each element is a single posterior sample.
Each matrix has four columns: \code{Column 1} identifies a node in \code{phy}. 
\code{Column 2} identifies the beginning time of the branch or segment 
of the branch that subtends the node in \code{Column 1}. \code{Column 3} identifies the 
ending time of the branch or segment of the branch that subtends the node in \code{Column 1}. 
\code{Column 4} identifies the index of the event that obtains along the branch or segment 
of the branch that subtends the node in \code{Column 1}.
}
\item{tipStates}{
A list of integer vectors. Each element is a single posterior sample.
For each tip the index of the event that obtains along the branch subtending the tip.
Tips are ordered increasing here and elsewhere.
}
\item{tipLambda}{
A list of numeric vectors. Each element is a single posterior sample.
For each tip the rate of speciation or trait evolution at the end of 
the terminal branch subtending that tip.
}
\item{tipMu}{
A list of numeric vectors. Each element is a single posterior sample.
For each tip the rate of extinction at the end of that terminal branch
subtending that tip. Meaningless if working with BAMM trait results.
}
\item{meanTipLambda}{
For each tip the mean of the marginal posterior density of the rate of speciation
or trait evolution at the end of each terminal branch subtending that tip.
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
the descendants of any node. e.g. if \code{downseq[2]=89} and \code{lastvist[89]=15}, the descendants 
of node 89 can be found by \code{downseq[2:15]}.
}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dan Rabosky, Mike Grundler
}
\note{
Currently the function does not check for duplicate tip labels in \code{phy}, which may cause the function to choke.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary.bammdata}}, \code{\link{plot.bammdata}}, \code{\link{dtRates}}.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
