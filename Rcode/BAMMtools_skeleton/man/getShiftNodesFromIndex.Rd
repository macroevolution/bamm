\name{getShiftNodesFromIndex}
\alias{getShiftNodesFromIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getShiftNodesFromIndex(ephy, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
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
function (ephy, index) 
{
    if (index > length(ephy$eventData)) {
        cat("Error. Attempting to access non-existent element from BAMM-data object\n")
        cat("You have << ", length(ephy$eventData), " >>> samples in your BAMM-data object\n")
        cat("Value of index must be no greater than the number of samples\n\n")
        stop()
    }
    root <- length(ephy$tip.label) + 1
    nodes <- ephy$eventData[[index]]$node
    nodes <- nodes[nodes != root]
    return(nodes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
