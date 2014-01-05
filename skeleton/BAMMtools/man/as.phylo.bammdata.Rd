\name{as.phylo.bammdata}
\alias{as.phylo.bammdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
as.phylo.bammdata(ephy)
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
    if (!"bammdata" \%in\% class(ephy)) {
        stop("Object ephy must be of class bammdata\n")
    }
    newphylo <- list()
    newphylo$edge <- ephy$edge
    newphylo$Nnode <- ephy$Nnode
    newphylo$tip.label <- ephy$tip.label
    newphylo$edge.length <- ephy$edge.length
    class(newphylo) <- "phylo"
    attributes(newphylo)$order = attributes(ephy)$order
    return(newphylo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
