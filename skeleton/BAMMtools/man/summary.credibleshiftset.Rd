\name{summary.credibleshiftset}
\alias{summary.credibleshiftset}
\alias{print.credibleshiftset}
\title{
Summary of credible set of shift configurations from a \code{BAMM} analysis.
}
\description{
Prints summary attributes of the \code{BAMM} credible set of shift configurations.
}
\usage{
\method{summary}{credibleshiftset}(object, ...)
\method{print}{credibleshiftset}(x, ...)
}
\arguments{
  \item{object,x}{
an object of class \code{credibleshiftset}.
}
\item{...}{
	additional arguments (unused).
}
}
\details{
Prints to console summary attributes of the XX\% credible set of shift configurations sampled using \code{BAMM}. Attributes printed include: the number of distinct configurations in the XX\% credible set and the posterior probability, cumulative probability, and number of rate shifts in the 9 most-probable shift configurations.
 
}
\value{
\code{summary.credibleshiftset} returns (invisibly) a dataframe with a number of rows equal to the number of shift configurations in the credible set and four columns:
\item{rank}{
The ranked index of each shift configuration (ranked by posterior probability)
}
\item{probability}{
the posterior probability of each shift configuration
}
\item{cumulative}{
the cumulative probability of each shift configuration
}
\item{N_shifts}{
the number of rate shifts in each shift configuration (can be zero).
}
}

\references{
\code{bamm-project.org}
}
\author{
Dan Rabosky
}

\seealso{
\code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammshifts}}, \code{\link{credibleShiftSet}}	
} 
 
 
\examples{
data(events.whales)
data(whales)
ed <- getEventData(whales, events.whales)
 
data(prior.whales)
bp <- getBranchShiftPriors(whales, prior.whales)
cset <- credibleShiftSet(ed, threshold = bp)
summary(cset)

}
 
\keyword{ models }
