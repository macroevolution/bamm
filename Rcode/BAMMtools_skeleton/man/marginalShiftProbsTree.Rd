\name{marginalShiftProbsTree}
\alias{marginalShiftProbsTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
marginalShiftProbsTree(ephy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ephy) 
{
    if (!"bamm-data" \%in\% class(ephy)) {
        stop("Object ephy must be of class bamm-data\n")
    }
    shiftvec <- numeric(length(ephy$edge.length))
    rootnode <- length(ephy$tip.label) + 1
    for (i in 1:length(ephy$eventData)) {
        hasShift <- ephy$edge[, 2] \%in\% ephy$eventData[[i]]$node
        shiftvec[hasShift] <- shiftvec[hasShift] + rep(1, sum(hasShift))
    }
    shiftvec <- shiftvec/length(ephy$eventData)
    newphy <- as.phylo.bammdata(ephy)
    newphy$edge.length <- shiftvec
    return(newphy)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
