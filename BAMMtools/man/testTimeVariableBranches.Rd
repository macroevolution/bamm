\name{testTimeVariableBranches}
\alias{testTimeVariableBranches}

\title{
Evaluate evidence for temporal rate variation across tree
}
\description{
For each branch in a phylogenetic tree, evaluates the evidence (posterior probability or Bayes factor) that macroevolutionary rates have varied through time 
}
\usage{
testTimeVariableBranches(ephy, prior_tv = 0.5, return.type = "posterior")
}

\arguments{
  \item{ephy}{
an object of class \code{bammdata}.
}
  \item{prior_tv}{
The prior probability that rate shifts lead to a new time-varying rate process (versus a time-constant process)
}
  \item{return.type}{
either \code{"posterior"} or \code{"bayesfactor"}, depending on which form of evidence you would like.
}
}
\details{
In \code{BAMM 2.0}, rate shifts on trees can lead to time-varying or constant-rate diversification processes. In other words, the model will incorporate temporal variation in rates only if there is sufficient evidence in the data to favor it. The function \code{testTimeVariableBranches} enables the user to extract the evidence in favor of time-varying rates on any branch of a phylogenetic tree from a \code{bammdata} object. 

The function returns a copy of the original phylogenetic tree, but where branch lengths have been replaced by either the posterior probability (\code{return.type = "posterior"}) or the Bayes factor evidence (\code{return.type = "bayesfactor"}) that the macroevolutionary rate regime governing each branch is time-variable. Consider a particular branch X on a phylogenetic tree. If the length of this branch is 0.97 and \code{return.type = "posterior"}, this implies that branch X was governed by a time-varying rate dynamic in 97\% of all samples in the posterior. Alternatively, only 3\% of samples specified a constant rate dynamic on this branch. 

The function also provides an alternative measure of support if \code{return.type = "bayesfactor"}. In this case, the Bayes factor evidence for temporal rate variation is computed for each branch. We simply imagine that diversification rates on each branch can be explained by one of two models: either rates vary through time, or they do not. In the above example (branch X), the Bayes factor would be computed as follows, letting \emph{Prob_timevar} and \emph{Prior_timevar} be the posterior and prior probabilities that a particular branch is governed by a time-varying rate process:

\emph{( Prob_timevar) / (1 - Prob_timevar)} * \emph{ (1 - prior_timevar) / (prior_timevar) }

The Bayes factor is not particularly useful under uniform prior odds (e.g., \code{prior_tv = 0.5}), since this simply reduces to the ratio of posterior probabilities. Note that the prior must correspond to whatever you used to analyze your data in \code{BAMM}. By default, time-variable and time-constant processes are assumed to have equal prior odds.

This function can be used several ways, but this function allows the user to quickly evaluate which portions of a phylogenetic tree have "significant" evidence for rate variation through time (see examples below). 

}

\author{
Dan Rabosky
}
\value{
An object of class \code{phylo}, but where branch lengths are replaced with the desired evidence (posterior probability or Bayes factor) that each branch is governed by a time-varying rate dynamic.
}


\examples{
# Load whale data:
data(whales)
data(events.whales)
ed <- getEventData(whales, events.whales, burnin=0.1, nsamples=200)

# compute the posterior probability of 
# time-varying rates on each branch
tree.pp <- testTimeVariableBranches(ed)

# Plot tree, but color all branches where the posterior 
# probability of time-varying rates exceeds 95\%:

colvec <- rep("black", nrow(whales$edge))
colvec[tree.pp$edge.length >= 0.95] <- 'red'

plot.phylo(whales, edge.color=colvec, cex=0.5)

# now, compute Bayes factors for each branch:

tree.bf <- testTimeVariableBranches(ed, return.type = "bayesfactor")

# now, assume that our prior was heavily stacked in favor
# of a time-constant process:
tree.bf2 <- testTimeVariableBranches(ed, prior_tv = 0.1, return.type = "bayesfactor")

# Plotting the branch-specific Bayes factors against each other:

plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0, 260), ylim=c(0, 260))
points(tree.bf2$edge.length, tree.bf$edge.length, pch=21, bg='red', cex=1.5)
axis(1)
axis(2, las=1)
mtext(side=1, text="Bayes factor: prior_tv = 0.1", line=3, cex=1.5)
mtext(side = 2, text = "Bayes factor: uniform prior odds", line=3, cex=1.5)

# and you can see that if your prior favors CONSTANT RATE dynamics
#  you will obtain much stronger Bayes factor support for time varying rates
# IF the evidence is present in your data to support time variation.
# To be clear, the Bayes factors in this example were computed from the 
#  same posterior probabilities: it is only the prior odds that differed.



}

