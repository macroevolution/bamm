\name{getSequenceForwardTraversal}
\alias{getSequenceForwardTraversal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSequenceForwardTraversal(phy, node)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
}
  \item{node}{
%%     ~~Describe \code{node} here~~
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
function (phy, node) 
{
    if (node <= length(phy$tip.label)) {
    }
    else {
        dset <- phy$edge[, 2][phy$edge[, 1] == node]
        phy$downseq <- c(phy$downseq, dset[1])
        phy <- getSequenceForwardTraversal(phy, dset[1])
        phy$downseq <- c(phy$downseq, dset[2])
        phy <- getSequenceForwardTraversal(phy, dset[2])
    }
    phy$lastvisit[node] <- phy$downseq[length(phy$downseq)]
    return(phy)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
